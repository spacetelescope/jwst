import logging
import re

from . import gwcs_drizzle
from . import resample_utils
from .. import datamodels
from ..lib.basic_utils import bytes2human
from ..model_blender import blendmeta

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["OutputTooLargeError", "ResampleData"]


class OutputTooLargeError(RuntimeError):
    """Raised when the output is too large for in-memory instantiation"""


class ResampleData:
    """
    This is the controlling routine for the resampling process.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model, such as pixfrac,
         weight type, exposure time (if relevant), and kernel, and merges
         them with any user-provided values.
      2. Creates output WCS based on input images and define mapping function
         between all input arrays and the output array.
      3. Updates output data model with output arrays from drizzle, including
         a record of metadata from all input models.
    """

    def __init__(self, input_models, output=None, single=False, blendheaders=True,
                 pixfrac=1.0, kernel="square", fillval="INDEF", weight_type="ivm",
                 good_bits=0, pscale_ratio=1.0, **kwargs):
        """
        Parameters
        ----------
        input_models : list of objects
            list of data models, one for each input image

        output : str
            filename for output

        kwargs : dict
            Other parameters
        """
        self.input_models = input_models

        self.pscale_ratio = pscale_ratio
        self.single = single
        self.blendheaders = blendheaders
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = weight_type
        self.good_bits = good_bits

        if output is None:
            output = input_models.meta.resample.output
        self.output_filename = output

        # Define output WCS based on all inputs, including a reference WCS
        self.output_wcs = resample_utils.make_output_wcs(self.input_models,
                                                         pscale_ratio=self.pscale_ratio)
        log.debug('Output mosaic size: {}'.format(self.output_wcs.data_size))
        can_allocate, required_memory = datamodels.util.check_memory_allocation(
            self.output_wcs.data_size, kwargs['allowed_memory'], datamodels.ImageModel
        )
        if not can_allocate:
            raise OutputTooLargeError(
                f'Combined ImageModel size {self.output_wcs.data_size} '
                f'requires {bytes2human(required_memory)}. '
                f'Model cannot be instantiated.'
            )
        self.blank_output = datamodels.ImageModel(self.output_wcs.data_size)

        # update meta data and wcs
        self.blank_output.update(input_models[0], only='PRIMARY')
        self.blank_output.update(input_models[0], only='SCI')
        self.blank_output.meta.wcs = self.output_wcs

        self.output_models = datamodels.ModelContainer()

    # TODO: Following method is not used, and never has been.  Find out if needed
    # def update_driz_outputs(self):
    #     """ Define output arrays for use with drizzle operations.
    #     """
    #     numchips = len(self.input_models)
    #     numplanes = (numchips // 32) + 1

    #     # Replace CONTEXT array with full set of planes needed for all inputs
    #     outcon = np.zeros((numplanes, self.output_wcs.data_size[0],
    #                        self.output_wcs.data_size[1]), dtype=np.int32)
    #     self.blank_output.con = outcon

    def do_drizzle(self):
        """Perform drizzling operation on input images's to create a new output
        """
        if self.single:
            return self.resample_many_to_many()
        else:
            return self.resample_many_to_one()

    def blend_output_metadata(self, output_model):
        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = output_model.meta.filename

        log.info('Blending metadata for {}'.format(output_file))
        blendmeta.blendmodels(output_model, inputs=self.input_models,
                              output=output_file)

    def resample_many_to_many(self):
        """Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure, i.e. map NRCA5 and
        NRCB5 onto the same output image, as they image different areas of the
        sky.

        Used for outlier detection
        """
        for exposure in self.input_models.models_grouped:
            output_model = self.blank_output.copy()
            output_model.update(exposure, only="PRIMARY")

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model, pixfrac=self.pixfrac,
                                            kernel=self.kernel, fillval=self.fillval)
            for img in exposure:
                # TODO: should weight_type=None here?
                inwht = resample_utils.build_driz_weight(img, weight_type=self.weight_type,
                                                         good_bits=self.good_bits)
                # apply sky subtraction
                blevel = img.meta.background.level
                if not img.meta.background.subtracted and blevel is not None:
                    data = img.data - blevel
                else:
                    data = img.data

                driz.add_image(data, img.meta.wcs, inwht=inwht)

            self.output_models.append(output_model)

        return self.output_models

    def resample_many_to_one(self):
        """Resample and coadd many inputs to a single output.

        Used for stage 3 resampling
        """
        output_model = self.blank_output.copy()
        output_model.meta.filename = self.output_filename
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pointings = len(self.input_models.group_names)

        if self.blendheaders:
            self.blend_output_metadata(output_model)

        # Initialize the output with the wcs
        driz = gwcs_drizzle.GWCSDrizzle(output_model, pixfrac=self.pixfrac,
                                        kernel=self.kernel, fillval=self.fillval)

        for img in self.input_models:
            inwht = resample_utils.build_driz_weight(img,
                                                     weight_type=self.weight_type,
                                                     good_bits=self.good_bits)
            # apply sky subtraction
            blevel = img.meta.background.level
            if not img.meta.background.subtracted and blevel is not None:
                data = img.data - blevel
            else:
                data = img.data

            driz.add_image(data, img.meta.wcs, inwht=inwht)

        self.update_fits_wcs(output_model)
        self.update_exposure_times(output_model)
        self.output_models.append(output_model)

        return self.output_models

    def update_exposure_times(self, output_model):
        """Modify exposure time metadata in-place"""
        total_exposure_time = 0.
        exposure_times = {'start': [], 'end': []}
        for exposure in self.input_models.models_grouped:
            total_exposure_time += exposure[0].meta.exposure.exposure_time
            exposure_times['start'].append(exposure[0].meta.exposure.start_time)
            exposure_times['end'].append(exposure[0].meta.exposure.end_time)

        # Update some basic exposure time values based on output_model
        output_model.meta.exposure.exposure_time = total_exposure_time
        output_model.meta.exposure.start_time = min(exposure_times['start'])
        output_model.meta.exposure.end_time = max(exposure_times['end'])
        output_model.meta.resample.product_exposure_time = total_exposure_time

    def update_fits_wcs(self, model):
        """
        Update FITS WCS keywords of the resampled image.
        """
        # Delete any SIP-related keywords first
        pattern = r"^(cd[12]_[12]|[ab]p?_\d_\d|[ab]p?_order)$"
        regex = re.compile(pattern)

        keys = list(model.meta.wcsinfo.instance.keys())
        for key in keys:
            if regex.match(key):
                del model.meta.wcsinfo.instance[key]

        # Write new PC-matrix-based WCS based on GWCS model
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
        model.meta.wcsinfo.ctype1 = "RA---TAN"
        model.meta.wcsinfo.ctype2 = "DEC--TAN"
