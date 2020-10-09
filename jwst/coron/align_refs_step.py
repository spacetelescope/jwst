""" Replace bad pixels and align psf image with target image."""

from jwst.datamodels.dqflags import interpret_bit_flags

from ..stpipe import Step
from .. import datamodels
from . import imageregistration
from . median_replace_img import median_replace_img

__all__ = ["AlignRefsStep"]


class AlignRefsStep(Step):

    """
    AlignRefsStep: Align coronagraphic PSF images
    with science target images.
    """

    spec = """
        median_box_length = integer(default=3,min=0) # box size for the median filter
        bad_bits = string(default="DO_NOT_USE") # the DQ bit values of bad pixels
    """

    reference_file_types = ['psfmask']

    def process(self, target, psf):

        # Open the input science target model
        with datamodels.open(target) as target_model:

            # Get the name of the psf mask reference file to use
            self.mask_name = self.get_reference_file(target_model, 'psfmask')
            self.log.info('Using PSFMASK reference file %s', self.mask_name)

            # Check for a valid reference file
            if self.mask_name == 'N/A':
                self.log.warning('No PSFMASK reference file found')
                self.log.warning('Align_refs step will be skipped')
                return None

            # Open the psf mask reference file
            mask_model = datamodels.ImageModel(self.mask_name)

            # Open the input psf images
            psf_model = datamodels.open(psf)

            # Retrieve the box size for the filter
            box_size = self.median_box_length

            # Get the bit value of bad pixels. A value of 0 treats all pixels as good.
            bad_bitvalue = self.bad_bits
            bad_bitvalue = interpret_bit_flags(bad_bitvalue)
            if bad_bitvalue is None:
                bad_bitvalue = 0

            # Replace bad pixels in the psf images
            psf_model = median_replace_img(psf_model, box_size, bad_bitvalue)

            # Replace bad pixels in the target images
            target_model = median_replace_img(target_model, box_size, bad_bitvalue)

            # Call the alignment routine
            result = imageregistration.align_models(target_model, psf_model,
                                                    mask_model)
            result.meta.cal_step.align_psfs = 'COMPLETE'

            mask_model.close()
            psf_model.close()

        return result
