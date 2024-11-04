import abc
import logging
import copy
import json
import math
import numpy as np

from typing import Union, Tuple, NamedTuple, List
from astropy.modeling import polynomial
from scipy.interpolate import CubicSpline
from scipy import interpolate
from scipy import ndimage
from astropy.io import fits
from gwcs import WCS

from stdatamodels.properties import ObjectNode
from stdatamodels.jwst import datamodels
from stdatamodels import DataModel
from stdatamodels.jwst.datamodels import dqflags, SlitModel, SpecModel,SpecwcsModel, MiriLrsPsfModel
from jwst.datamodels import SourceModelContainer

from ..assign_wcs.util import wcs_bbox_from_shape
from ..lib import pipe_utils
from ..lib.wcs_utils import get_wavelengths
from . import extract1d
from . import ifu
from . import spec_wcs
from . import utils   # these are routines taken from extract.py and pulled out of the class definition. 

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""


def open_specwcs(specwcs_ref_name: str, exp_type: str):
    """Open the specwcs reference file.

    Parameters
    ----------
    specwcs_ref_name : str
        The name of the specwcs reference file. This file contains information of the trace location.
        For the MIRI LRS FIED-SlIT it is a fits file containing the x,y center of the trace. 
    ext_type : str
        The exposure type of the data

    Returns
    -------
        Currently only works on MIRI LRS-FIXEDSLIT exposures. 
        Return the center of the trace in x and y for a given wavelength

    """
    
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        # use fits to read file (the datamodel does not have all that is needed) 
        ref = fits.open(specwcs_ref_name)

        with ref:
            lrsdata = np.array([d for d in ref[1].data])
            # Get the zero point from the reference data.
            # The zero_point is X, Y  (which should be COLUMN, ROW)
            # These are 1-indexed in CDP-7 (i.e., SIAF convention) so must be converted to 0-indexed
            # for lrs_fixedslit
            zero_point = ref[0].header['imx'] - 1, ref[0].header['imy'] - 1
    
        # In the lrsdata reference table, X_center,Y_center, wavelength  relative to zero_point

        xcen = lrsdata[:, 0]
        ycen = lrsdata[:, 1]
        wavetab = lrsdata[:, 2]
        trace = xcen + zero_point[0]
        wave_trace = ycen + zero_point[1]
    return trace, wave_trace, wavetab
    

def open_psf(psf_refname: str, exp_type: str):
    """Open the PSF reference file.

    Parameters
    ----------
    psf_ref_name : str
        The name of the psf reference file. 
    ext_type : str
        The exposure type of the data

    Returns
    -------
        Currenly only works on MIRI LRS-FIXEDSLIT exposures. 
        Returns the EPSF model.

    """
    
    psf_model = None
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        psf_model = MiriLrsPsfModel(psf_refname)
        # The information we read in from PSF file is:
        # center_col: psf_model.meta.psf.center_col
        # super sample factor: psf_model.meta.psf.subpix)
        # psf : psf_model.data (2d)
        # wavelength of PSF planes: psf_model.wave 
    return psf_model 

    

# TODO JEM  THIS CLASS IS NOT BEING USED YET. JUST HOOKS FOR LATER
class PSFExtractData():
    """ Class for PSF extraction for each source"""
    
    def __init__(self,
                 input_model: DataModel,
                 slit: Union[DataModel, None] = None,
                 dispaxis: int = HORIZONTAL):

        """
        Parameters
        ----------
        input_model : data model
            The input science data.

        slit : an input slit, or None if not used
            For MultiSlit, `slit` is one slit from
            a list of slits in the input.  For other types of data, `slit`
            will not be used.

        dispaxis : int
            Dispersion direction:  1 is horizontal, 2 is vertical.
        """
  
        self.exp_type = input_model.meta.exposure.type        
        self.dispaxis = dispaxis
        self.wcs = None

        if slit is None:
            if hasattr(input_model.meta, 'wcs'):
                self.wcs = input_model.meta.wcs
        elif hasattr(slit, 'meta') and hasattr(slit.meta, 'wcs'):
            self.wcs = slit.meta.wcs

        if self.wcs is None:
            log.warning("WCS function not found in input.")


def psf_extraction(input_model, psf_ref_name, specwcs_ref_name, save_spatial_profile):
    """Open the PSF reference file.

    Parameters
    ----------
    input_model : data model
        This can be either the input science file or one SlitModel out of
        a list of slits.
    psf_ref_name : str
        PSF reference filename
    specwcs_ref_name : str
        Reference file containing information on trace of spectra

    save_spatial_profile : bool
        Save the spatial profile. 
    Returns
    -------
        Currenly only works on MIRI LRS-FIXEDSLIT exposures. 
        spatial profile

    """
        
    exp_type = input_model.meta.exposure.type
    trace, wave_trace, wavetab = open_specwcs(specwcs_ref_name, exp_type)
    psf_model = open_psf(psf_ref_name, exp_type)
    # Input data 
    # 1. Must bePoint Source (for now PSF extraction is only working with Point Source data)
    # 2. Can have multiple sources in one input_model
    # 3. Keep track of mode to help decide how to find the sources and use WCS
    
    # We want to organize the data by source. The input can be a source ModelContainer or a single model.
    # The type of the individual models can be: ImageModel (LRS), SlitModel (NIRSpec Slit), MultiSlitModel (MOS)
    # of IFUCubeModel
    # For now store all the data  as a list of input models. We will loop over each input mode, find the sources
    # determine the spatial profile, and call extract1d. 
    
    input_models = []
    if isinstance(input_model, SourceModelContainer):
        input_models = input_model
    else:
        input_models.append(input_model)
        
    for model in input_models:

        # currently for LRS-FIXEDSLIT data it is assumes there is 1 source and the source is located
        # at the dither_position
        
        if exp_type == 'MIR_LRS-FIXEDSLIT': # assume we can use the WCS to find source location
            dispaxis = 2
            wcs  = model.meta.wcs

            # assuming 1 source
            # Set up the profiles that are to be determined
            nsource = 1
            profile_2d = []
            profile_bg = []  # Currently not sure if this is needed (at least for MIRI LRS data) 
                             # Calwebb_spec2 subtracts the other nod position automatically 

            # Only one source determine the location using the WCS
            middle, locn_w, locn_x = utils.locn_from_wcs(
                dispaxis,
                wcs, 
                model,
                None, None, None)

            # Perform fit of reference trace and corresponding wavelength
            # The wavelength for the reference trace does not exactly line up exactly with the data
            
            cs = CubicSpline(wavetab, trace)
            cen_shift = cs(locn_w)
            shift = locn_x - cen_shift

            log.info ('Location of source in observed spectrum %s at wavelength %s', locn_x, locn_w)

            log.info('For this wavelength: %s  the reference trace the location is %s ', locn_w, cen_shift)
            log.info('Shift to apply to ref trace %s', shift)
            
            # For LRS-Fix data  we want to cut out the slit from the full array
            slit = []
            bbox = wcs.bounding_box
            x0, x1 = bbox[0]
            y0, y1 = bbox[1]
            cutout = model.data[int(np.ceil(y0)):int(np.ceil(y1)), int(np.round(x0)):int(np.round(x1))]
            slit.append(cutout)
            
            # adjust the trace to the slit region 
            trace_cutout = trace - bbox[0][0]
            wtrace_cutout = wave_trace - bbox[1][0]
            trace_shift = trace_cutout + shift
            psf_wave = psf_model.wave
            
            # trace_shift - for each wavelength in the PSF this is the shift in x to apply to the PSF image to shift it
            #                 to fall on the source.
            # wavetab : is the wavelengh cooresponding the the trace. This wavelength may not match exactly to the the PSF.
            # 
            # determine what the shifts per row are for the the wavelengths given by the model PSF 

            psf_subpix = psf_model.meta.psf.subpix
            
            PSFinterp = interpolate.interp1d(wavetab, trace_shift,fill_value="extrapolate")
            psf_shift = PSFinterp(psf_wave)
            psf_shift =  psf_model.meta.psf.center_col - (psf_shift * psf_subpix)

            # if the PSF is an ePSF
            _x, _y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))

            # Convert 1D psf to 2 d with second dimension = 1
            psf2 = psf_shift[:, np.newaxis]
            _x = _x*psf_subpix + psf2

            sprofile = ndimage.map_coordinates(psf_model.data, [_y, _x], order=1)
            profile_2d.append(sprofile)

            # TODO JEM - this might be just for testing. If we keep writing spatial profile
            # we will need a datamodel module for it  rather than using FITS write
            # Also need to change how we form the name if we keep writing spatial profile
            if save_spatial_profile:
                file_out = model.meta.filename.replace('.fits', '_spatial_profile.fits')
                hdu = fits.PrimaryHDU()
                hdu.header['SUBPIX'] = psf_model.meta.psf.subpix
                hdu.header['CENTCOL'] = psf_model.meta.psf.center_col
                hdu.header['CENTWAVE'] = locn_w
                hdu.header['CENTLOC'] = locn_x
                
                hdu.header['REFSHFT'] = float(cen_shift)
                hdu.header['TRCSHFT'] = shift
                hdu1 = fits.ImageHDU(sprofile, name='SPROFILE')
                hdu_final = fits.HDUList([hdu, hdu1])

                print(file_out) 
                hdu_final.writeto(file_out, overwrite=True)

            # set up the data to be passed to extract1d
            # TODO JEM - setup data - gathers the data that extract1d will need.
            # we do not have to do this in seperate module and we need to know what
            # Tim B routine needs - so this will need to be updated. 
            result = utils.setup_data(model)
            data, dq, var_poisson, var_rnoise, weights = result
            






