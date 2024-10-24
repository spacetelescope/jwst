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
from stdatamodels.jwst.datamodels import dqflags, SlitModel, SpecModel,SpecwcsModel, PsfModel
from jwst.datamodels import SourceModelContainer

from ..assign_wcs.util import wcs_bbox_from_shape
from ..lib import pipe_utils
from ..lib.wcs_utils import get_wavelengths
from . import extract1d
from . import ifu
from . import spec_wcs
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)



HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""

# This is intended to be larger than any possible distance (in pixels) between the target and any point in the image;
# used by locn_from_wcs().
HUGE_DIST = 1.e20

def open_specwcs(specwcs_ref_name: str, exp_type):
    print('in open specwcs', exp_type)
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        # use fits to read file (the datamodel does not have all that is needed) 
        ref = fits.open(specwcs_ref_name)
        print(' Number of rows in the table', ref[1].data.shape)

        with ref:
            lrsdata = np.array([d for d in ref[1].data])
            # Get the zero point from the reference data.
            # The zero_point is X, Y  (which should be COLUMN, ROW)
            # These are 1-indexed in CDP-7 (i.e., SIAF convention) so must be converted to 0-indexed
            # for lrs_fixedslit
            zero_point = ref[0].header['imx'] - 1, ref[0].header['imy'] - 1
            print('zero point', zero_point) 
    
        # In the lrsdata reference table, X_center,Y_center, wavelength  relative to zero_point
        # x0,y0(ul) x1,y1 (ur) x2,y2(lr) x3,y3(ll) define corners of the box within which the distortion
        # and wavelength calibration was derived
        xcen = lrsdata[:, 0]
        ycen = lrsdata[:, 1]
        wavetab = lrsdata[:, 2]
        trace = xcen + zero_point[0]
        wave_trace = ycen + zero_point[1]
    return trace, wave_trace, wavetab
    

def open_psf(psf_refname: str, exp_type):

    psf_model = None
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        psf_model = PsfModel(psf_refname)
        print(psf_model.meta.psf.center_col)
        print(psf_model.meta.psf.subpix)
        print(psf_model.data)
        print(psf_model.wave) 
    return psf_model 

    

class PSFExtractData():
    """ Class for PSF extraction for each source"""
    
    def __init__(self,
                 input_model: DataModel,
                 slit: Union[DataModel, None] = None,
                 dispaxis: int = HORIZONTAL,
                 use_source_posn: Union[bool, None] = None):


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



            
            
def psf_extraction(input_model, psf_ref_name, specwcs_ref_name):

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

        # currently for LRS-FIXEDSLIT data it is assume there is 1 source and the source is located
        # at the dither_position
        
        if exp_type == 'MIR_LRS-FIXEDSLIT': # assume we can use the WCS to find source location
            dispaxis = 2
            wcs  = model.meta.wcs

            # assuming 1 source
            # Set up the profiles that are to be determined
            nsource = 1
            profile_2d = []
            profile_bg = []  # Currently not sure if this is needed
                             # Calwebb_spec2 subtracts the other nod position automatically 

            # Only one source 
            middle, locn_w, locn_x = utils.locn_from_wcs(
                dispaxis,
                wcs, 
                model,
                None, None, None)

            # Match the  reference trace and corresponding wavelength to the data
            # The wavelength for the reference trace does not exactly line up exactly with the data
            
            cs = CubicSpline(wavetab, trace)
            cen_shift = cs(locn_w)
            shift = locn_x - cen_shift

            print('Location of source in observed spectrum', locn_x, ' at wavelength', locn_w)
            print('In the reference trace the location is ', cen_shift)
            print('Shift to apply to ref trace', shift)
            
            # For LRS-Fix data  we want to cut out the slit from the full array
            slit = []
            bbox = wcs.bounding_box
            x0, x1 = bbox[0]
            y0, y1 = bbox[1]
            cutout = model.data[int(np.ceil(y0)):int(np.ceil(y1)), int(np.round(x0)):int(np.round(x1))]
            slit.append(cutout)
            
            # adjust the trace to for slit region 
            trace_cutout = trace - bbox[0][0]
            wtrace_cutout = wave_trace - bbox[1][0]
            trace_shift = trace_cutout + xshift
   
            # xtrace_shift - for each wavelength in the PSF this is the shift in x to apply to the PSF image to shift it
            #                 to fall on the source.
            # wavetab : is the wavelengh cooresponding the the trace. This wavelength may not match exactly to the the PSF.wave

            # determine what the shifts per row are for the the wavelengths given by the model PSF 

            PSFinterp = interpolate.interp1d(wavetab, xtrace_shift,fill_value="extrapolate")
            psf_shift = PSFinterp(psf_model.wave)
            psf_shift =  psf_model.meta.center_col - (psf_shift * psf_subpix)

            # if the PSF is an ePSF
            _x, _y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
            _x = _x*psf_subpix + psf_shift[:, np.newaxis]
            sprofile = ndimage.map_coordinates(psf, [_y, _x], order=1)
            profile_2d.append(sprofile)
            
            print(sprofile)

            # set up the data to be passed to extract1d
            result = util.setup_data(model)
            data, dq, var_poisson, var_rnoise, weights = result
            






