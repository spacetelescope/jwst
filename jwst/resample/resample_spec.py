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

from .. import datamodels
from . import gwcs_drizzle
from . import resample_utils


CRBIT = np.uint32(datamodels.dqflags.pixel['JUMP_DET'])

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

    def __init__(self, input_models, output=None, **pars):
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

        self.drizpars = pars
        self.pscale_ratio = 1.
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
        Create a spatial/spectral WCS output frame

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

        if refmodel is None:
            refmodel = self.input_models[0]

        refwcs = refmodel.meta.wcs
        bb = refwcs.bounding_box
        grid = wcstools.grid_from_bounding_box(bb)
        ra, dec, lam = np.array(refwcs(*grid))
        lon = np.nanmean(ra)
        lat = np.nanmean(dec)
        tan = Pix2Sky_TAN()
        native2celestial = RotateNative2Celestial(lon, lat, 180)
        undist2sky = tan | native2celestial
        # Filter out RuntimeWarnings due to computed NaNs in the WCS
        warnings.simplefilter("ignore")
        x_tan, y_tan = undist2sky.inverse(ra, dec)
        warnings.resetwarnings()

        spectral_axis = find_dispersion_axis(refmodel)
        spatial_axis = spectral_axis ^ 1

        # Compute the wavelength array, trimming NaNs from the ends
        wavelength_array = np.nanmedian(lam, axis=spectral_axis)
        wavelength_array = wavelength_array[~np.isnan(wavelength_array)]

        # Compute RA and Dec up the slit (spatial direction) at the center
        # of the dispersion.  Use spectral_axis to determine slicing dimension
        lam_center_index = int((bb[spectral_axis][1] - bb[spectral_axis][0]) / 2)
        if not spectral_axis:
            x_tan_array = x_tan.T[lam_center_index]
            y_tan_array = y_tan.T[lam_center_index]
        else:
            x_tan_array = x_tan[lam_center_index]
            y_tan_array = y_tan[lam_center_index]
        x_tan_array = x_tan_array[~np.isnan(x_tan_array)]
        y_tan_array = y_tan_array[~np.isnan(y_tan_array)]

        fitter = LinearLSQFitter()
        fit_model = Linear1D()
        pix_to_ra = fitter(fit_model, np.arange(x_tan_array.shape[0]), x_tan_array)
        pix_to_dec = fitter(fit_model, np.arange(y_tan_array.shape[0]), y_tan_array)

        # Tabular interpolation model, pixels -> lambda
        pix_to_wavelength = Tabular1D(lookup_table=wavelength_array,
            bounds_error=False, fill_value=None, name='pix2wavelength')

        # Tabular models need an inverse explicitly defined.
        # If the wavelength array is decending instead of ascending, both
        # points and lookup_table need to be reversed in the inverse transform
        # for scipy.interpolate to work properly
        points = wavelength_array
        lookup_table = np.arange(wavelength_array.shape[0])
        if not np.all(np.diff(wavelength_array) > 0):
            points = points[::-1]
            lookup_table = lookup_table[::-1]
        pix_to_wavelength.inverse = Tabular1D(points=points,
            lookup_table=lookup_table,
            bounds_error=False, fill_value=None, name='wavelength2pix')

        # For the input mapping, duplicate the spatial coordinate
        mapping = Mapping((spatial_axis, spatial_axis, spectral_axis))

        # Sometimes the slit is perpendicular to the RA or Dec axis.
        # For example, if the slit is perpendicular to RA, that means
        # the slope of pix_to_ra will be nearly zero, so make sure
        # mapping.inverse uses pix_to_dec.inverse.  The auto definition
        # of mapping.inverse is to use the 2nd spatial coordinate, i.e. Dec.
        if np.isclose(pix_to_dec.slope, 0, atol=1e-8):
            mapping_tuple = (0, 1)
            # Account for vertical or horizontal dispersion on detector
            if spatial_axis:
                mapping.inverse = Mapping(mapping_tuple[::-1])
            else:
                mapping.inverse = Mapping(mapping_tuple)

        # The final transform
        transform = mapping | (pix_to_ra & pix_to_dec | undist2sky) & pix_to_wavelength

        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
            reference_frame=coord.ICRS())
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,),
            unit=(u.micron,), axes_names=('wavelength',))
        world = cf.CompositeFrame([sky, spec], name='world')

        pipeline = [(det, transform),
                    (world, None)]

        output_wcs = WCS(pipeline)

        # compute the output array size in WCS axes order, i.e. (x, y)
        output_array_size = [0, 0]
        output_array_size[spectral_axis] = len(wavelength_array)
        output_array_size[spatial_axis] = len(x_tan_array)

        # turn the size into a numpy shape in (y, x) order
        self.data_size = tuple(output_array_size[::-1])
        bounding_box = resample_utils.wcs_bbox_from_shape(self.data_size)
        output_wcs.bounding_box = bounding_box

        return output_wcs

    def do_drizzle(self, xmin=0, xmax=0, ymin=0, ymax=0, **pars):
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

            bb = resample_utils.wcs_bbox_from_shape(output_model.data.shape)
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
                    weight_type=self.drizpars['weight_type'],
                    good_bits=self.drizpars['good_bits'])

                if hasattr(img, 'name'):
                    log.info('Resampling slit {} {}'.format(img.name, self.data_size))
                else:
                    log.info('Resampling slit {}'.format(self.data_size))

                in_wcs = img.meta.wcs
                driz.add_image(img.data, in_wcs, inwht=inwht,
                               expin=img.meta.exposure.exposure_time,
                               pscale_ratio=self.pscale_ratio,
                               xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

            # Update some basic exposure time values based on all the inputs
            output_model.meta.exposure.exposure_time = texptime
            output_model.meta.exposure.start_time = min(exposure_times['start'])
            output_model.meta.exposure.end_time = max(exposure_times['end'])
            output_model.meta.resample.product_exposure_time = texptime
            output_model.meta.resample.weight_type = self.drizpars['weight_type']
            output_model.meta.resample.pointings = pointings

            # Update mutlislit slit info on the output_model
            for attr in ['name', 'xstart', 'xsize', 'ystart', 'ysize',
                         'slitlet_id', 'source_id', 'source_name', 'source_alias',
                         'stellarity', 'source_type', 'source_xpos', 'source_ypos',
                         'dispersion_direction', 'shutter_state']:
                try:
                    val = getattr(img, attr)
                except AttributeError:
                    pass
                else:
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
