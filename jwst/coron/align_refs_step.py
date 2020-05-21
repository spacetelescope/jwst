""" Smooth and align psf image with target image."""
from ..stpipe import Step
from .. import datamodels
from . import imageregistration
from . smooth_psf import smooth_psf

__all__ = ["AlignRefsStep"]


class AlignRefsStep(Step):

    """
    AlignRefsStep: Align coronagraphic PSF images
    with science target images.
    """

    spec = """
        smoothing_box_length = integer(default=4,min=0) # Smoothing kernel
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

            # Retrieve the lenth of the smoothing box for the filter
            try:
                box_size = self.smoothing_length
            except AttributeError:
                box_size = 4
                self.log.warning('Smoothing length not found, set to default')

            # Smooth the psf images
            psf_model = smooth_psf(psf_model, box_size)

            # Call the alignment routine
            result = imageregistration.align_models(target_model, psf_model,
                                                    mask_model)
            result.meta.cal_step.align_psfs = 'COMPLETE'

            mask_model.close()
            psf_model.close()

        return result
