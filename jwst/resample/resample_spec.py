from __future__ import (division, print_function, unicode_literals, 
    absolute_import)

import time
from collections import OrderedDict
import numpy as np
from scipy import interpolate

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.modeling.models import (Shift, Scale, Mapping, Rotation2D,
    Pix2Sky_TAN, RotateNative2Celestial)

from .. import datamodels
from ..assign_wcs import util
from . import gwcs_drizzle
from . import bitmask
from . import resample_utils
from gwcs import wcstools, WCS
from gwcs.utils import _compute_lon_pole

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResampleSpecData(object):
    """
    This is the controlling routine for the resampling process.
    It loads and sets the various input data and parameters needed by
    the drizzle function and then calls the C-based cdriz.tdriz function
    to do the actual resampling.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model, such as pixfrac,
         weight type, exposure time (if relevant), and kernel, and merges
         them with any user-provided values.
      2. Creates output WCS based on input images and define mapping function
         between all input arrays and the output array.
      3. Initializes all output arrays, including WHT and CTX arrays.
      4. Passes all information for each input chip to drizzle function.
      5. Updates output data model with output arrays from drizzle, including
         (eventually) a record of metadata from all input models.
    """
    drizpars = {'single': False,
                'kernel': 'square',
                'pixfrac': 1.0,
                'good_bits': 4,
                'fillval': 'INDEF',
                'wht_type': 'exptime'}

    def __init__(self, input_models, output=None, ref_filename=None, **pars):
        """
        Parameters
        ----------
        input_models : list of objects
            list of data models, one for each input image

        output : str
            filename for output
        """
        self.input_models = input_models
        if output is None:
            output = input_models.meta.resample.output
        self.output_filename = output
        self.ref_filename = ref_filename

        self.output_spatial_scale = None
        self.output_spectral_scale = None
        self.output_wcs = None

        # If user specifies use of drizpars ref file (default for pipeline use)
        # update input parameters with default values from ref file
        if self.ref_filename is not None:
            self.get_drizpars()
        self.drizpars.update(pars)

        # Define output WCS based on all inputs, including a reference WCS
        wcslist = [m.meta.wcs for m in self.input_models]
        self.build_output_wcs()
        self.build_size_from_domain()
        self.blank_output = datamodels.DrizProductModel(self.data_size)
        self.blank_output.assign_wcs(self.output_wcs)

        # Default to defining output models metadata as
        # a copy of the first input_model's metadata
        ### TO DO:
        ###    replace this with a call to a generalized version of fitsblender
        ###
        self.blank_output.update(self.input_models[0])
        self.output_models = datamodels.ModelContainer()


    def get_drizpars(self):
        """ Extract drizzle parameters from reference file
        """
        # start by interpreting input data models to define selection criteria
        num_groups = len(self.input_models.group_names)
        input_dm = self.input_models[0]
        filtname = input_dm.meta.instrument.filter

        # Create a data model for the reference file
        ref_model = datamodels.DrizParsModel(self.ref_filename)
        # look for row that applies to this set of input data models
        # NOTE:
        # This logic could be replaced by a method added to the DrizParsModel 
        # object to select the correct row based on a set of selection params
        row = None
        drizpars = ref_model.drizpars_table

        filter_match = False # flag to support wild-card rows in drizpars table
        for n, filt, num in zip(range(1, drizpars.numimages.shape[0] + 1), 
            drizpars.filter, drizpars.numimages):
            # only remember this row if no exact match has already been made for
            # the filter. This allows the wild-card row to be anywhere in the
            # table; since it may be placed at beginning or end of table.

            if filt == b"ANY" and not filter_match and num_groups >= num:
                row = n
            # always go for an exact match if present, though...
            if filtname == filt and num_groups >= num:
                row = n
                filter_match = True

        # With presence of wild-card rows, code should never trigger this logic
        if row is None:
            txt = "No row found in {0} that matches input data."
            log.error(txt.format(self.ref_filename))
            raise ValueError

        # read in values from that row for each parameter
        for kw in self.drizpars:
            if kw in drizpars.names:
                self.drizpars[kw] = ref_model['drizpars_table.{0}'.format(kw)][row]


    def create_output_metadata(self):
        """ Create new output metadata based on blending all input metadata
        """
        # TODO: Modify API for fitsblender
        pass


    def build_output_wcs(self, refwcs=None):
        """
        Create a simple output wcs covering footprint of the input datamodels
        """
        # TODO: generalize this for more than one input datamodel
        # TODO: generalize this for imaging modes with distorted wcs
        input_model = self.input_models[0]
        if refwcs == None:
            refwcs = input_model.meta.wcs

        # Generate grid of sky coordinates for area within domain
        det = x, y = wcstools.grid_from_domain(refwcs.domain)
        sky = ra, dec, lam = refwcs(*det)
        domain_xsize = refwcs.domain[0]['upper'] - refwcs.domain[0]['lower']
        domain_ysize = refwcs.domain[1]['upper'] - refwcs.domain[1]['lower']
        x_center, y_center = int(domain_xsize / 2), int(domain_ysize / 2)

        # Compute slit angular size, slit center sky coords
        xpos = []
        sz = 3
        for row in lam:
            if np.isnan(row[x_center]):
                xpos.append(np.nan)
            else:
                f = interpolate.interp1d(row[x_center-sz+1:x_center+sz],
                    x[y_center, x_center-sz+1:x_center+sz])
                xpos.append(f(lam[y_center, x_center]))
        x_arg = np.array(xpos)[~np.isnan(lam[:,x_center])]
        y_arg = y[~np.isnan(lam[:,x_center]),x_center]
        # slit_coords, spect0 = refwcs(x_arg, y_arg, output='numericals_plus')
        slit_ra, slit_dec, slit_spec_ref = refwcs(x_arg, y_arg)
        slit_coords = SkyCoord(ra=slit_ra, dec=slit_dec, unit=u.deg)
        pix_num = np.flipud(np.arange(len(slit_ra)))
        # pix_num = np.arange(len(slit_ra))
        interpol_ra = interpolate.interp1d(pix_num, slit_ra)
        interpol_dec = interpolate.interp1d(pix_num, slit_dec)
        slit_center_pix = len(slit_spec_ref) / 2. - 1
        log.debug('Slit center pix: {0}'.format(slit_center_pix))
        slit_center_sky = SkyCoord(ra=interpol_ra(slit_center_pix),
            dec=interpol_dec(slit_center_pix), unit=u.deg)
        log.debug('Slit center: {0}'.format(slit_center_sky))
        log.debug('Fiducial: {0}'.format(resample_utils.compute_spec_fiducial([refwcs])))
        angular_slit_size = np.abs(slit_coords[0].separation(slit_coords[-1]))
        log.debug('Slit angular size: {0}'.format(angular_slit_size.arcsec))
        dra, ddec = slit_coords[0].spherical_offsets_to(slit_coords[-1])
        offset_up_slit = (dra.to(u.arcsec), ddec.to(u.arcsec))
        log.debug('Offset up the slit: {0}'.format(offset_up_slit))

        # Compute spatial and spectral scales
        xposn = np.array(xpos)[~np.isnan(xpos)]
        dx = xposn[-1] - xposn[0]
        slit_npix = np.sqrt(dx**2 + np.array(len(xposn) - 1)**2)
        spatial_scale = angular_slit_size / slit_npix
        log.debug('Spatial scale: {0}'.format(spatial_scale.arcsec))
        spectral_scale = lam[y_center, x_center] - lam[y_center, x_center - 1]

        # Compute slit angle relative (clockwise) to y axis
        slit_rot_angle = (np.arcsin(dx / slit_npix) * u.radian).to(u.degree)
        log.debug('Slit rotation angle: {0}'.format(slit_rot_angle))

        # Compute transform for output frame
        roll_ref = input_model.meta.wcsinfo.roll_ref * u.deg
        min_lam = np.nanmin(lam)
        offset = Shift(0.) & Shift(slit_center_pix * -1.)
        # TODO: double-check the signs on the following rotation angles
        # rot = Rotation2D(roll_ref)
        rot = Rotation2D(roll_ref + slit_rot_angle)
        scale = Scale(spatial_scale) & Scale(spatial_scale)
        tan = Pix2Sky_TAN()
        lon_pole = _compute_lon_pole(slit_center_sky, tan)
        skyrot = RotateNative2Celestial(slit_center_sky.ra, slit_center_sky.dec,
            lon_pole)
        spatial_trans = offset | rot | scale | tan | skyrot
        spectral_trans = Scale(spectral_scale) | Shift(min_lam)
        mapping = Mapping((1, 1, 0),)
        mapping.inverse = Mapping((2, 1))
        transform = mapping | spatial_trans & spectral_trans
        transform.outputs = ('ra', 'dec', 'lamda')

        # Build the output wcs
        input_frame = refwcs.input_frame
        output_frame = refwcs.output_frame
        wnew = WCS(output_frame=output_frame, forward_transform=transform)

        # Build the domain in the output frame wcs object
        domain_grid = wnew.backward_transform(*sky)
        domain = []
        for axis in input_frame.axes_order:
            axis_min = np.nanmin(domain_grid[axis])
            axis_max = np.nanmax(domain_grid[axis]) + 1
            domain.append({'lower': axis_min, 'upper': axis_max,
                'includes_lower': True, 'includes_upper': False})
        log.debug('Domain: {0} {1}'.format(domain[1]['lower'], domain[1]['upper']))
        wnew.domain = domain

        # Update class properties
        self.output_spatial_scale = spatial_scale
        self.output_spectral_scale = spectral_scale
        self.output_wcs = wnew


    def build_size_from_domain(self, refwcs=None):
        """ Compute the size of the output frame based on the domain
        """
        if refwcs == None:
            refwcs = self.output_wcs
        size = []
        for axis in refwcs.domain:
            delta = axis['upper'] - axis['lower']
            size.append(int(delta + 0.5))
        self.data_size = tuple(reversed(size))
        log.debug('DATA SIZE {0}'.format(self.data_size))


    def do_drizzle(self, **pars):
        """ Perform drizzling operation on input images's to create a new output
        """
        # Set up information about what outputs we need to create: single or final
        # Key: value from metadata for output/observation name
        # Value: full filename for output file
        driz_outputs = OrderedDict()

        # Look for input configuration parameter telling the code to run
        # in single-drizzle mode (mosaic all detectors in a single observation?)
        if self.drizpars['single']:
            driz_outputs = ['{0}_resamp.fits'.format(g) for g in self.input_models.group_names]
            model_groups = self.input_models.models_grouped
            group_exptime = []
            for group in model_groups:
                group_exptime.append(group[0].meta.exposure.exposure_time)
        else:
            final_output = self.input_models.meta.resample.output # get global name
            driz_outputs = [final_output]
            model_groups = [self.input_models]

            total_exposure_time = 0.0
            for group in self.input_models.models_grouped:
                total_exposure_time += group[0].meta.exposure.exposure_time
            group_exptime = [total_exposure_time]

        pointings = len(self.input_models.group_names)
        # Now, generate each output for all input_models
        for obs_product, group, texptime in zip(driz_outputs, model_groups, group_exptime):
            output_model = self.blank_output.copy()
            output_model.meta.wcs = self.output_wcs
            # # TODO: do this properly in wcs_from_spec_footprints()
            # output_model.meta.wcs.domain = self.output_wcs.domain
            # # Instead we do a cludge below to get the domain to not be neg.
            output_model.meta.wcs.domain = resample_utils.create_domain(
                self.output_wcs, output_model.data.shape)
            output_model.meta.filename = obs_product

            output_model.meta.asn.pool_name = self.input_models.meta.pool_name
            output_model.meta.asn.table_name = self.input_models.meta.table_name

            exposure_times = {'start': [], 'end': []}

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model,
                                single=self.drizpars['single'],
                                pixfrac=self.drizpars['pixfrac'],
                                kernel=self.drizpars['kernel'],
                                fillval=self.drizpars['fillval'])

            for n, img in enumerate(group):
                exposure_times['start'].append(img.meta.exposure.start_time)
                exposure_times['end'].append(img.meta.exposure.end_time)

                # outwcs_pscale = output_model.meta.wcsinfo.cdelt3
                # wcslin_pscale = img.meta.wcsinfo.cdelt3
                outwcs_pixel_scale = self.output_spatial_scale
                inwcs_pixel_scale = self.output_spatial_scale
                pscale_ratio = outwcs_pixel_scale / inwcs_pixel_scale

                inwht = build_driz_weight(img, wht_type=self.drizpars['wht_type'],
                                    good_bits=self.drizpars['good_bits'])
                log.info('Resampling slit {0} {1}'.format(img.name,
                    self.data_size))
                driz.add_image(img.data, img.meta.wcs, inwht=inwht,
                        expin=img.meta.exposure.exposure_time,
                        pscale_ratio=pscale_ratio)

            # Update some basic exposure time values based on all the inputs
            output_model.meta.exposure.exposure_time = texptime
            output_model.meta.exposure.start_time = min(exposure_times['start'])
            output_model.meta.exposure.end_time = max(exposure_times['end'])
            output_model.meta.resample.product_exposure_time = texptime
            output_model.meta.resample.product_data_extname = driz.sciext
            output_model.meta.resample.product_context_extname = driz.conext
            output_model.meta.resample.product_weight_extname = driz.whtext
            output_model.meta.resample.drizzle_fill_value = str(driz.fillval)
            output_model.meta.resample.drizzle_pixel_fraction = driz.pixfrac
            output_model.meta.resample.drizzle_kernel = driz.kernel
            output_model.meta.resample.drizzle_output_units = driz.out_units
            output_model.meta.resample.drizzle_weight_scale = driz.wt_scl
            output_model.meta.resample.resample_bits = self.drizpars['good_bits']
            output_model.meta.resample.weight_type = self.drizpars['wht_type']
            output_model.meta.resample.pointings = pointings

            # Update mutlislit slit info on the output_model
            for attr in ['name', 'xstart', 'xsize', 'ystart', 'ysize',
                'slitlet_id', 'source_id', 'source_name', 'source_alias',
                'catalog_id', 'stellarity', 'source_type', 'source_xpos',
                'source_ypos', 'nshutters', 'relsens']:
                if hasattr(img, attr):
                    setattr(output_model, attr, getattr(img, attr))

            self.output_models.append(output_model)

        #self.output_models.save(None)  # DEBUG: Remove for production


def build_mask(dqarr, bitvalue):
    """ Builds a bit-mask from an input DQ array and a bitvalue flag
    """

    bitvalue = bitmask.interpret_bits_value(bitvalue)

    if bitvalue is None:
        return (np.ones(dqarr.shape, dtype=np.uint8))
    return np.logical_not(np.bitwise_and(dqarr, ~bitvalue)).astype(np.uint8)


def build_driz_weight(model, wht_type=None, good_bits=None):
    """ Create input weighting image based on user inputs
    """
    if good_bits is not None and good_bits < 0:
        good_bits = None
    dqmask = build_mask(model.dq, good_bits)
    exptime = model.meta.exposure.exposure_time

    if wht_type.lower()[:3] == 'err':
        # Multiply the scaled ERR file by the input mask in place.
        inwht = (exptime / model.err)**2 * dqmask
    #elif wht_type == 'IVM':
    #    _inwht = img.buildIVMmask(chip._chip,dqarr,pix_ratio)
    elif wht_type.lower()[:3] == 'exp':# or wht_type is None, as used for single=True
        inwht = exptime * dqmask
    else:
        # Create an identity input weight map
        inwht = np.ones(model.data.shape, dtype=model.data.dtype)
    return inwht
