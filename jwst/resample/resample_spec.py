import logging
from collections import OrderedDict
import warnings
import numpy as np

from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import (Mapping, Tabular1D, Linear1D,
                                     Pix2Sky_TAN, RotateNative2Celestial)
from astropy.modeling.fitting import LinearLSQFitter
from gwcs import wcstools, WCS
from gwcs import coordinate_frames as cf
from ..cube_build.cube_build_wcs_util import  wrap_ra

from .. import datamodels
from . import gwcs_drizzle
from . import resample_utils
from ..model_blender import blendmeta


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

    def __init__(self, input_models, output=None, single=False,
        blendheaders=False, pixfrac=1.0, kernel="square", fillval=0,
        weight_type="exptime", good_bits=0, pscale_ratio=1.0, **kwargs):
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

        self.pscale_ratio = pscale_ratio
        self.single = single
        self.blendheaders = blendheaders
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = weight_type
        self.good_bits = good_bits

        self.blank_output = None

        # Define output WCS based on all inputs, including a reference WCS
        self.output_wcs = self.build_interpolated_output_wcs()
        self.blank_output = datamodels.SlitModel(self.data_size)
        self.blank_output.update(self.input_models[0], only="PRIMARY")
        self.blank_output.update(self.input_models[0], only="SCI")
        self.blank_output.meta.wcs = self.output_wcs
        self.output_models = datamodels.ModelContainer()

    def build_interpolated_output_wcs(self, refmodel=None):
        """
        Create a spatial/spectral WCS output frame using all the input models

        Creates output frame by linearly fitting RA, Dec along the slit and
        producing a lookup table to interpolate wavelengths in the dispersion
        direction.

        Parameters
        ----------
        refmodel : `~jwst.datamodels.DataModel`
            The reference input image from which the fiducial WCS is created.
            If not specified, the first image in self.input_models is used.

        Returns
        -------
        output_wcs : `~gwcs.WCS` object
            A gwcs WCS object defining the output frame WCS
        """

        # for each input model convert slit x,y to ra,dec,lam
        # use first input model to set spatial scale
        # use center of appended ra and dec arrays to set up
        # center of final ra,dec
        # append all ra,dec, wavelength array for each slit
        # use first model to initialize wavelenth array
        # append wavelengths that fall outside the endpoint of
        # of wavelength array when looping over additional data

        all_wavelength = []
        all_ra_slit = []
        all_dec_slit = []

        for im, model in enumerate(self.input_models):
            wcs = model.meta.wcs
            bb = wcs.bounding_box
            grid = wcstools.grid_from_bounding_box(bb)
            ra, dec, lam = np.array(wcs(*grid))
            spectral_axis = find_dispersion_axis(model)
            spatial_axis = spectral_axis ^ 1

            # Compute the wavelength array, trimming NaNs from the ends
            # In many cases, a whole slice is NaNs, so ignore those warnings
            warnings.simplefilter("ignore")
            wavelength_array = np.nanmedian(lam, axis=spectral_axis)
            warnings.resetwarnings()
            wavelength_array = wavelength_array[~np.isnan(wavelength_array)]

            # We need to estimate the spatial sampling to use for the output WCS.
            # Tt is assumed the spatial sampling is the same for all the input
            # models. So we can use the first input model to set the spatial
            # sampling.

            # Steps to do this for first input model:
            # 1. find the middle of the spectrum in wavelength
            # 2. Pull out the ra and dec at the center of the slit.
            # 3. Find the mean ra,dec and the center of the slit this will
            #    represent the tangent point
            # 4. Convert ra,dec -> tangent plane projection: x_tan,y_tan
            # 5. using x_tan, y_tan perform a linear fit to find spatial sampling
            # first input model sets intializes wavelength array and defines
            # the spatial scale of the output wcs
            if im == 0:
                for iw in wavelength_array:
                    all_wavelength.append(iw)

                lam_center_index = int((bb[spectral_axis][1] -
                                            bb[spectral_axis][0]) / 2)
                if spatial_axis == 0:
                    ra_center = ra[lam_center_index,:]
                    dec_center = dec[lam_center_index,:]
                else:
                    ra_center = ra[:,lam_center_index]
                    dec_center = dec[:,lam_center_index]
                # find the ra and dec for this slit using center of slit
                ra_center_pt = np.nanmean(ra_center)
                dec_center_pt = np.nanmean(dec_center)

                if resample_utils.is_sky_like(model.meta.wcs.output_frame):
                    # convert ra and dec to tangent projection
                    tan = Pix2Sky_TAN()
                    native2celestial = RotateNative2Celestial(ra_center_pt, dec_center_pt, 180)
                    undist2sky1 = tan | native2celestial
                    # Filter out RuntimeWarnings due to computed NaNs in the WCS
                    warnings.simplefilter("ignore")
                    # at this center of slit find x,y tangent projection - x_tan, y_tan
                    x_tan, y_tan = undist2sky1.inverse(ra, dec)
                    warnings.resetwarnings()
                else:
                    # for non sky-like output frames, no need to do tangent plane projections
                    # but we still use the same variables
                    x_tan, y_tan = ra, dec

                # pull out data from center
                if spectral_axis == 0:
                    x_tan_array = x_tan.T[lam_center_index]
                    y_tan_array = y_tan.T[lam_center_index]
                else:
                    x_tan_array = x_tan[lam_center_index]
                    y_tan_array = y_tan[lam_center_index]

                x_tan_array = x_tan_array[~np.isnan(x_tan_array)]
                y_tan_array = y_tan_array[~np.isnan(y_tan_array)]

                # estimate the spatial sampling
                fitter = LinearLSQFitter()
                fit_model = Linear1D()
                xstop = x_tan_array.shape[0] / self.pscale_ratio
                xstep = 1 / self.pscale_ratio
                ystop = y_tan_array.shape[0] / self.pscale_ratio
                ystep = 1 / self.pscale_ratio
                pix_to_xtan = fitter(fit_model, np.arange(0, xstop, xstep), x_tan_array)
                pix_to_ytan = fitter(fit_model, np.arange(0, ystop, ystep), y_tan_array)

            # append all ra and dec values to use later to find min and max
            # ra and dec
            ra_use = ra.flatten()
            ra_use = ra_use[~np.isnan(ra_use)]
            dec_use = dec.flatten()
            dec_use = dec_use[~np.isnan(dec_use)]
            all_ra_slit.append(ra_use)
            all_dec_slit.append(dec_use)

            # now check wavelength array to see if we need to add to it
            this_minw = np.min(wavelength_array)
            this_maxw = np.max(wavelength_array)
            all_minw = np.min(all_wavelength)
            all_maxw = np.max(all_wavelength)

            if this_minw < all_minw:
                addpts = wavelength_array[wavelength_array < all_minw]
                for ip in range(len(addpts)):
                    all_wavelength.append(addpts[ip])
            if this_maxw > all_maxw:
                addpts = wavelength_array[wavelength_array > all_maxw]
                for ip in range(len(addpts)):
                    all_wavelength.append(addpts[ip])

        # done looping over set of models

        all_ra = np.hstack(all_ra_slit)
        all_dec = np.hstack(all_dec_slit)
        all_wave = np.hstack(all_wavelength)
        all_wave = all_wave[~np.isnan(all_wave)]
        all_wave = np.sort(all_wave, axis=None)
        # Tabular interpolation model, pixels -> lambda
        wavelength_array = np.unique(all_wave)
        # Check if the data is MIRI LRS FIXED Slit. If it is then
        # the wavelength array needs to be flipped so that the resampled
        # dispersion direction matches the disperion direction on the detector.
        if self.input_models[0].meta.exposure.type == 'MIR_LRS-FIXEDSLIT':
            wavelength_array = np.flip(wavelength_array, axis=None)

        step = 1 / self.pscale_ratio
        stop = wavelength_array.shape[0] / self.pscale_ratio
        points = np.arange(0, stop, step)
        pix_to_wavelength = Tabular1D(points=points,
                                      lookup_table=wavelength_array,
                                      bounds_error=False, fill_value=None,
                                      name='pix2wavelength')

        # Tabular models need an inverse explicitly defined.
        # If the wavelength array is decending instead of ascending, both
        # points and lookup_table need to be reversed in the inverse transform
        # for scipy.interpolate to work properly
        points = wavelength_array
        lookup_table = np.arange(0, stop, step)

        if not np.all(np.diff(wavelength_array) > 0):
            points = points[::-1]
            lookup_table = lookup_table[::-1]
        pix_to_wavelength.inverse = Tabular1D(points=points,
                                              lookup_table=lookup_table,
                                              bounds_error=False, fill_value=None,
                                              name='wavelength2pix')

        # For the input mapping, duplicate the spatial coordinate
        mapping = Mapping((spatial_axis, spatial_axis, spectral_axis))

        # Sometimes the slit is perpendicular to the RA or Dec axis.
        # For example, if the slit is perpendicular to RA, that means
        # the slope of pix_to_xtan will be nearly zero, so make sure
        # mapping.inverse uses pix_to_ytan.inverse.  The auto definition
        # of mapping.inverse is to use the 2nd spatial coordinate, i.e. Dec.

        if np.isclose(pix_to_ytan.slope, 0, atol=1e-8):
            mapping_tuple = (0, 1)
            # Account for vertical or horizontal dispersion on detector
            if spatial_axis:
                mapping.inverse = Mapping(mapping_tuple[::-1])
            else:
                mapping.inverse = Mapping(mapping_tuple)

        # The final transform
        # redefine the ra, dec center tangent point to include all data

        # check if all_ra crosses 0 degress - this makes it hard to
        # define the min and max ra correctly
        all_ra = wrap_ra(all_ra)
        ra_min = np.amin(all_ra)
        ra_max = np.amax(all_ra)

        ra_center_final  = (ra_max + ra_min)/2.0
        dec_min = np.amin(all_dec)
        dec_max = np.amax(all_dec)
        dec_center_final  = (dec_max + dec_min)/2.0
        tan = Pix2Sky_TAN()
        if len(self.input_models) == 1: # single model use ra_center_pt to be consistent
                                        # with how resample was done before
            ra_center_final = ra_center_pt
            dec_center_final = dec_center_pt

        if resample_utils.is_sky_like(model.meta.wcs.output_frame):
            native2celestial = RotateNative2Celestial(ra_center_final, dec_center_final, 180)
            undist2sky = tan | native2celestial
            # find the spatial size of the output - same in x,y
            x_tan_all, _ = undist2sky.inverse(all_ra, all_dec)
        else:
            x_tan_all, _ = all_ra, all_dec
        x_min = np.amin(x_tan_all)
        x_max = np.amax(x_tan_all)
        x_size = int(np.ceil((x_max - x_min)/np.absolute(pix_to_xtan.slope)))

        # single model use size of x_tan_array
        # to be consistent with method before
        if len(self.input_models) == 1:
            x_size = len(x_tan_array)

        # define the output wcs
        if resample_utils.is_sky_like(model.meta.wcs.output_frame):
            transform = mapping | (pix_to_xtan & pix_to_ytan | undist2sky) & pix_to_wavelength
        else:
            transform = mapping | (pix_to_xtan & pix_to_ytan) & pix_to_wavelength

        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        if resample_utils.is_sky_like(model.meta.wcs.output_frame):
            sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                                    reference_frame=coord.ICRS())
        else:
            sky = cf.Frame2D(name=f'resampled_{model.meta.wcs.output_frame.name}', axes_order=(0, 1))
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,),
                                unit=(u.micron,), axes_names=('wavelength',))
        world = cf.CompositeFrame([sky, spec], name='world')

        pipeline = [(det, transform),
                    (world, None)]

        output_wcs = WCS(pipeline)

        # import ipdb; ipdb.set_trace()

        # compute the output array size in WCS axes order, i.e. (x, y)
        output_array_size = [0, 0]
        output_array_size[spectral_axis] = int(np.ceil(len(wavelength_array) / self.pscale_ratio))
        output_array_size[spatial_axis] = int(np.ceil(x_size / self.pscale_ratio))
        # turn the size into a numpy shape in (y, x) order
        self.data_size = tuple(output_array_size[::-1])
        bounding_box = resample_utils.wcs_bbox_from_shape(self.data_size)
        output_wcs.bounding_box = bounding_box

        return output_wcs

    def blend_output_metadata(self, output_model):
        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = output_model.meta.filename

        log.info(f'Blending metadata for {output_file}')
        blendmeta.blendmodels(output_model, inputs=self.input_models, output=output_file)

    def do_drizzle(self, xmin=0, xmax=0, ymin=0, ymax=0, **pars):
        """ Perform drizzling operation on input images's to create a new output
        """
        # Set up information about what outputs we need to create: single or final
        # Key: value from metadata for output/observation name
        # Value: full filename for output file
        driz_outputs = OrderedDict()

        # Look for input configuration parameter telling the code to run
        # in single-drizzle mode (mosaic all detectors in a single observation?)
        if self.single:
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

            bb = resample_utils.wcs_bbox_from_shape(output_model.data.shape)
            output_model.meta.wcs.bounding_box = bb
            output_model.meta.filename = obs_product

            if self.blendheaders:
                self.blend_output_metadata(output_model)

            exposure_times = {'start': [], 'end': []}

            outwcs = output_model.meta.wcs

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model,
                                outwcs=outwcs,
                                single=self.single,
                                pixfrac=self.pixfrac,
                                kernel=self.kernel,
                                fillval=self.fillval)

            for n, img in enumerate(group):
                exposure_times['start'].append(img.meta.exposure.start_time)
                exposure_times['end'].append(img.meta.exposure.end_time)
                inwht = resample_utils.build_driz_weight(img,
                    weight_type=self.weight_type,
                    good_bits=self.good_bits)

                if hasattr(img, 'name'):
                    log.info('Resampling slit {} {}'.format(img.name, self.data_size))
                else:
                    log.info('Resampling slit {}'.format(self.data_size))

                in_wcs = img.meta.wcs
                driz.add_image(img.data, in_wcs, inwht=inwht,
                               expin=img.meta.exposure.exposure_time,
                               xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

            # Update some basic exposure time values based on all the inputs
            output_model.meta.exposure.exposure_time = texptime
            output_model.meta.exposure.start_time = min(exposure_times['start'])
            output_model.meta.exposure.end_time = max(exposure_times['end'])
            output_model.meta.resample.product_exposure_time = texptime
            output_model.meta.resample.weight_type = self.weight_type
            output_model.meta.resample.pointings = pointings

            # Update slit info on the output_model. This is needed
            # because model.slit attributes are not in model.meta, so the
            # normal update() method doesn't work with them.
            for attr in ['name', 'xstart', 'xsize', 'ystart', 'ysize',
                         'slitlet_id', 'source_id', 'source_name', 'source_alias',
                         'stellarity', 'source_type', 'source_xpos', 'source_ypos',
                         'dispersion_direction', 'shutter_state']:
                try:
                    val = getattr(img, attr)
                except AttributeError:
                    pass
                else:
                    if val is not None:
                        setattr(output_model, attr, val)

            self.output_models.append(output_model)

        return self.output_models


def find_dispersion_axis(refmodel):
    """
    Find the dispersion axis (0-indexed) of the given 2D wavelength array
    """
    dispaxis = refmodel.meta.wcsinfo.dispersion_direction
    # Change from 1 --> X and 2 --> Y to 0 --> X and 1 --> Y.
    return dispaxis - 1
