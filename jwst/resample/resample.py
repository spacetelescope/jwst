import logging
from collections import OrderedDict
import numpy as np

from .. import datamodels

from . import gwcs_drizzle
from . import resample_utils
from ..model_blender import blendmeta

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleData"]


class ResampleData:
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
        self.drizpars = pars
        if output is None:
            output = input_models.meta.resample.output
        self.output_filename = output

        # Define output WCS based on all inputs, including a reference WCS
        self.output_wcs = resample_utils.make_output_wcs(self.input_models)
        log.debug('Output mosaic size: {}'.format(self.output_wcs.data_size))
        self.blank_output = datamodels.ImageModel(self.output_wcs.data_size)

        # update meta data and wcs
        self.blank_output.update(input_models[0], only='PRIMARY')
        self.blank_output.update(input_models[0], only='SCI')
        self.blank_output.meta.wcs = self.output_wcs

        self.output_models = datamodels.ModelContainer()

    def update_driz_outputs(self):
        """ Define output arrays for use with drizzle operations.
        """
        numchips = len(self.input_models)
        numplanes = (numchips // 32) + 1

        # Replace CONTEXT array with full set of planes needed for all inputs
        outcon = np.zeros((numplanes, self.output_wcs.data_size[0],
                           self.output_wcs.data_size[1]), dtype=np.int32)
        self.blank_output.con = outcon

    def blend_output_metadata(self, output_model):
        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = output_model.meta.filename

        log.info('Blending metadata for {}'.format(output_file))
        blendmeta.blendmodels(output_model, inputs=self.input_models,
                              output=output_file)

    def do_drizzle(self):
        """ Perform drizzling operation on input images's to create a new output
        """
        # Set up information about what outputs we need to create: single or final
        # Key: value from metadata for output/observation name
        # Value: full filename for output file
        driz_outputs = OrderedDict()

        # Look for input configuration parameter telling the code to run
        # in single-drizzle mode (mosaic all detectors in a single observation)
        if self.drizpars['single']:
            driz_outputs = self.input_models.group_names
            exposures = self.input_models.models_grouped
            group_exptime = []
            for exposure in exposures:
                group_exptime.append(exposure[0].meta.exposure.exposure_time)
        else:
            driz_outputs = [self.output_filename]
            exposures = [self.input_models]

            total_exposure_time = 0.0
            for exposure in exposures:
                total_exposure_time += exposure[0].meta.exposure.exposure_time
            group_exptime = [total_exposure_time]
        pointings = len(self.input_models.group_names)

        for obs_product, exposure, texptime in zip(driz_outputs, exposures,
                                                   group_exptime):
            output_model = self.blank_output.copy()
            output_model.meta.filename = obs_product
            saved_model_type = output_model.meta.model_type

            if self.drizpars['blendheaders']:
                self.blend_output_metadata(output_model)
                output_model.meta.model_type = saved_model_type

            exposure_times = {'start': [], 'end': []}

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model,
                                            single=self.drizpars['single'],
                                            pixfrac=self.drizpars['pixfrac'],
                                            kernel=self.drizpars['kernel'],
                                            fillval=self.drizpars['fillval'])

            for n, img in enumerate(exposure):
                exposure_times['start'].append(img.meta.exposure.start_time)
                exposure_times['end'].append(img.meta.exposure.end_time)

                # apply sky subtraction
                blevel = img.meta.background.level
                if not img.meta.background.subtracted and blevel is not None:
                    img.data -= blevel

                outwcs_pscale = output_model.meta.wcsinfo.cdelt1
                wcslin_pscale = img.meta.wcsinfo.cdelt1

                inwht = resample_utils.build_driz_weight(img,
                    weight_type=self.drizpars['weight_type'],
                    good_bits=self.drizpars['good_bits'])
                driz.add_image(img.data, img.meta.wcs, inwht=inwht,
                        expin=img.meta.exposure.exposure_time,
                        pscale_ratio=outwcs_pscale / wcslin_pscale)

            # Update some basic exposure time values based on all the inputs
            output_model.meta.exposure.exposure_time = texptime
            output_model.meta.exposure.start_time = min(exposure_times['start'])
            output_model.meta.exposure.end_time = max(exposure_times['end'])
            output_model.meta.resample.product_exposure_time = texptime
            output_model.meta.resample.weight_type = self.drizpars['weight_type']
            output_model.meta.resample.pointings = pointings

            self.update_fits_wcs(output_model)

            self.output_models.append(output_model)

    def update_fits_wcs(self, model):
        """
        Update FITS WCS keywords of the resampled image.
        """
        transform = model.meta.wcs.forward_transform
        model.meta.wcsinfo.crpix1 = -transform[0].offset.value + 1
        model.meta.wcsinfo.crpix2 = -transform[1].offset.value + 1
        model.meta.wcsinfo.cdelt1 = transform[3].factor.value
        model.meta.wcsinfo.cdelt2 = transform[4].factor.value
        model.meta.wcsinfo.ra_ref = transform[6].lon.value
        model.meta.wcsinfo.dec_ref = transform[6].lat.value
        model.meta.wcsinfo.crval1 = model.meta.wcsinfo.ra_ref
        model.meta.wcsinfo.crval2 = model.meta.wcsinfo.dec_ref
        model.meta.wcsinfo.pc1_1 = transform[2].matrix.value[0][0]
        model.meta.wcsinfo.pc1_2 = transform[2].matrix.value[0][1]
        model.meta.wcsinfo.pc2_1 = transform[2].matrix.value[1][0]
        model.meta.wcsinfo.pc2_2 = transform[2].matrix.value[1][1]
