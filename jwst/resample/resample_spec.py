import functools
from collections import OrderedDict, namedtuple

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import (Shift, Scale, Mapping, Identity,
    Rotation2D, Pix2Sky_TAN, RotateNative2Celestial, Tabular1D, Const1D)
from gwcs import wcstools, WCS
from gwcs.utils import _compute_lon_pole
from gwcs import coordinate_frames as cf

from .. import datamodels
from . import gwcs_drizzle
from . import resample_utils

CRBIT = np.uint32(datamodels.dqflags.pixel['JUMP_DET'])

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResampleSpecData:
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
                'good_bits': CRBIT,
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
        self.ref_filename = ref_filename

        self.pscale_ratio = 1.
        self.blank_output = None

        # If user specifies use of drizpars ref file (default for pipeline use)
        # update input parameters with default values from ref file
        if self.ref_filename is not None:
            self.get_drizpars()
        self.drizpars.update(pars)

        # Define output WCS based on all inputs, including a reference WCS
        # wcslist = [m.meta.wcs for m in self.input_models]

        self.instrument_name = self.input_models[0].meta.instrument.name
        if self.instrument_name in ['NIRSPEC']:
            self.output_wcs = self.build_nirspec_output_wcs()
        elif self.instrument_name in ['MIRI']:
            self.output_wcs = self.build_miri_output_wcs()

        self.data_size = self.build_size_from_bounding_box()
        self.blank_output = datamodels.DrizProductModel(self.data_size)

        self.blank_output.update(datamodels.ImageModel(self.input_models[0]._instance))
        self.blank_output.meta.wcs = self.output_wcs
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

        # Flag to support wild-card rows in drizpars table
        filter_match = False
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


    def build_nirspec_output_wcs(self, refmodel=None):
        """
        Create a spatial/spectral WCS covering footprint of the input
        """
        if not refmodel:
            refmodel = self.input_models[0]
        refwcs = refmodel.meta.wcs
        bb = refwcs.bounding_box

        ref_det2slit = refwcs.get_transform('detector', 'slit_frame')
        ref_slit2world = refwcs.get_transform('slit_frame', 'world')

        grid = x, y = wcstools.grid_from_bounding_box(bb, step=(1, 1))
        Grid = namedtuple('Grid', refwcs.slit_frame.axes_names)
        grid_slit = Grid(*ref_det2slit(*grid))

        # Compute spatial transform from detector to slit
        ref_wavelength = np.nanmean(grid_slit.wavelength)

        # find the number of pixels sampled by a single shutter
        fid = np.array([[0., 0.], [-.5, .5], np.repeat(ref_wavelength, 2)])
        slit_extents = np.array(ref_det2slit.inverse(*fid)).T
        pix_per_shutter = np.linalg.norm(slit_extents[0] - slit_extents[1])

        # Get min and max of slit in pixel units
        ymin = np.nanmin(grid_slit.y_slit) * pix_per_shutter
        ymax = np.nanmax(grid_slit.y_slit) * pix_per_shutter
        slit_height_pix = int(abs(ymax - ymin + 0.5))

        # Compute grid of wavelengths and make tabular model w/inverse
        lookup_table = np.nanmean(grid_slit.wavelength, axis=0)
        wavelength_transform = Tabular1D(lookup_table=lookup_table,
            bounds_error=False, fill_value=np.nan)
        wavelength_transform.inverse = Tabular1D(points=lookup_table,
            lookup_table=np.arange(grid_slit.wavelength.shape[1]),
            bounds_error=False, fill_value=np.nan)

        # Define detector to slit transforms
        yslit_transform = Scale(-1 / pix_per_shutter) | Shift(ymax / pix_per_shutter)
        xslit_transform = Const1D(-0.)
        xslit_transform.inverse = Const1D(0.)

        # Construct the final transform
        coord_mapping = Mapping((0, 1, 0))
        coord_mapping.inverse = Mapping((2, 1))
        the_transform = xslit_transform & yslit_transform & wavelength_transform
        out_det2slit = coord_mapping | the_transform

        # Create coordinate frames
        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        slit_spatial = cf.Frame2D(name='slit_spatial', axes_order=(0, 1),
            unit=("", ""), axes_names=('x_slit', 'y_slit'))
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,),
            unit=(u.micron,), axes_names=('wavelength',))
        slit_frame = cf.CompositeFrame([slit_spatial, spec], name='slit_frame')
        sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
            reference_frame=coord.ICRS())
        world = cf.CompositeFrame([sky, spec], name='world')

        pipeline = [(det, out_det2slit),
                    (slit_frame, ref_slit2world),
                    (world, None)]

        output_wcs = WCS(pipeline)

        self.data_size = (slit_height_pix, len(lookup_table))
        bounding_box = resample_utils.bounding_box_from_shape(self.data_size)
        output_wcs.bounding_box = bounding_box
        return output_wcs


    def build_miri_output_wcs(self, refmodel=None):
        """
        Create a simple output wcs covering footprint of the input datamodels
        """
        # TODO: generalize this for more than one input datamodel
        # TODO: generalize this for imaging modes with distorted wcs
        if not refmodel:
            refmodel = self.input_models[0]
        refwcs = refmodel.meta.wcs

        x, y = wcstools.grid_from_bounding_box(refwcs.bounding_box, step=(1, 1),
            center=True)
        ra, dec, lam = refwcs(x.flatten(), y.flatten())
        # TODO: once astropy.modeling._Tabular is fixed, take out the
        # flatten() and reshape() code above and below
        ra = ra.reshape(x.shape)
        dec = dec.reshape(x.shape)
        lam = lam.reshape(x.shape)

        # Find rotation of the slit from y axis from the wcs forward transform
        # TODO: figure out if angle is necessary for MIRI.  See for discussion
        # https://github.com/STScI-JWST/jwst/pull/347
        rotation = [m for m in refwcs.forward_transform
                    if isinstance(m, Rotation2D)]
        if rotation:
            rot_slit = functools.reduce(lambda x, y: x | y, rotation)
            unrotate = rot_slit.inverse
            refwcs_minus_rot = refwcs.forward_transform | \
                unrotate & Identity(1)
            # Correct for this rotation in the wcs
            ra, dec, lam = refwcs_minus_rot(x.flatten(), y.flatten())
            ra = ra.reshape(x.shape)
            dec = dec.reshape(x.shape)
            lam = lam.reshape(x.shape)

        # Get the slit size at the center of the dispersion
        sky_coords = SkyCoord(ra=ra, dec=dec, unit=u.deg)
        slit_coords = sky_coords[int(sky_coords.shape[0] / 2)]
        slit_angular_size = slit_coords[0].separation(slit_coords[-1])
        log.debug('Slit angular size: {0}'.format(slit_angular_size.arcsec))

        # Compute slit center from bounding_box
        dx0 = refwcs.bounding_box[0][0]
        dx1 = refwcs.bounding_box[0][1]
        dy0 = refwcs.bounding_box[1][0]
        dy1 = refwcs.bounding_box[1][1]
        slit_center_pix = (dx1 - dx0) / 2
        dispersion_center_pix = (dy1 - dy0) / 2
        slit_center = refwcs_minus_rot(dx0 + slit_center_pix, dy0 +
            dispersion_center_pix)
        slit_center_sky = SkyCoord(ra=slit_center[0], dec=slit_center[1],
            unit=u.deg)
        log.debug('slit center: {0}'.format(slit_center))

        # Compute spatial and spectral scales
        spatial_scale = slit_angular_size / slit_coords.shape[0]
        log.debug('Spatial scale: {0}'.format(spatial_scale.arcsec))
        tcenter = int((dx1 - dx0) / 2)
        trace = lam[:, tcenter]
        trace = trace[~np.isnan(trace)]
        spectral_scale = np.abs((trace[-1] - trace[0]) / trace.shape[0])
        log.debug('spectral scale: {0}'.format(spectral_scale))

        # Compute transform for output frame
        log.debug('Slit center %s' % slit_center_pix)
        # TODO: double-check the signs on the following rotation angles
        roll_ref = refmodel.meta.wcsinfo.roll_ref * u.deg
        rot = Rotation2D(roll_ref)
        tan = Pix2Sky_TAN()
        lon_pole = _compute_lon_pole(slit_center_sky, tan)
        skyrot = RotateNative2Celestial(slit_center_sky.ra, slit_center_sky.dec,
            lon_pole)
        min_lam = np.nanmin(lam)
        mapping = Mapping((0, 0, 1))

        transform = Shift(-slit_center_pix) & Identity(1) | \
            Scale(spatial_scale) & Scale(spectral_scale) | \
            Identity(1) & Shift(min_lam) | mapping | \
            (rot | tan | skyrot) & Identity(1)

        transform.inputs = (x, y)
        transform.outputs = ('ra', 'dec', 'lamda')

        # Build the output wcs
        input_frame = refwcs.input_frame
        output_frame = refwcs.output_frame
        wnew = WCS(output_frame=output_frame, forward_transform=transform)

        # Build the bounding_box in the output frame wcs object
        bounding_box_grid = wnew.backward_transform(ra, dec, lam)

        bounding_box = []
        for axis in input_frame.axes_order:
            axis_min = np.nanmin(bounding_box_grid[axis])
            axis_max = np.nanmax(bounding_box_grid[axis])
            bounding_box.append((axis_min, axis_max))
        wnew.bounding_box = tuple(bounding_box)
        return wnew


    def build_size_from_bounding_box(self, refwcs=None):
        """ Compute the size of the output frame based on the bounding_box
        """
        if not refwcs:
            refwcs = self.output_wcs
        size = []
        for axis in refwcs.bounding_box:
            delta = axis[1] - axis[0]
            size.append(int(delta + 0.5))
        return tuple(reversed(size))


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
            final_output = self.input_models.meta.resample.output
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

            bb = resample_utils.bounding_box_from_shape(output_model.data.shape)
            output_model.meta.wcs.bounding_box = bb
            output_model.meta.filename = obs_product

            exposure_times = {'start': [], 'end': []}

            outwcs = output_model.meta.wcs

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model,
                                outwcs=outwcs,
                                single=self.drizpars['single'],
                                pixfrac=self.drizpars['pixfrac'],
                                kernel=self.drizpars['kernel'],
                                fillval=self.drizpars['fillval'])

            for n, img in enumerate(group):
                exposure_times['start'].append(img.meta.exposure.start_time)
                exposure_times['end'].append(img.meta.exposure.end_time)

                inwht = resample_utils.build_driz_weight(img,
                    wht_type=self.drizpars['wht_type'],
                    good_bits=self.drizpars['good_bits'])
                if hasattr(img, 'name'):
                    log.info('Resampling slit {} {}'.format(img.name, self.data_size))
                else:
                    log.info('Resampling slit {}'.format(self.data_size))

                in_wcs = img.meta.wcs
                driz.add_image(img.data, in_wcs, inwht=inwht,
                        expin=img.meta.exposure.exposure_time,
                        pscale_ratio=self.pscale_ratio)

            # Update some basic exposure time values based on all the inputs
            output_model.meta.exposure.exposure_time = texptime
            output_model.meta.exposure.start_time = min(exposure_times['start'])
            output_model.meta.exposure.end_time = max(exposure_times['end'])
            output_model.meta.resample.product_exposure_time = texptime
            output_model.meta.resample.weight_type = self.drizpars['wht_type']
            output_model.meta.resample.pointings = pointings

            # Update mutlislit slit info on the output_model
            try:
                for attr in ['name', 'xstart', 'xsize', 'ystart', 'ysize',
                        'slitlet_id', 'source_id', 'source_name', 'source_alias',
                        'stellarity', 'source_type', 'source_xpos', 'source_ypos',
                        'shutter_state', 'relsens']:
                    try:
                        val = getattr(img, attr)
                    except:
                        val = None
                    if val is not None:
                        setattr(output_model, attr, val)
            except:
                pass

            self.output_models.append(output_model)
