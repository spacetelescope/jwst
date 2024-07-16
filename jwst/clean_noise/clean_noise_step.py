from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import clean_noise

__all__ = ["CleanNoiseStep"]


class CleanNoiseStep(Step):
    """
    CleanNoiseStep: This step performs 1/f noise correction ("cleaning").

    Input data is expected to be a ramp file, in between jump and
    ramp fitting steps.

    Algorithms implemented are:
        - `fft`: Background noise is fit in frequency space.
           Implementation is based on the NSClean algorithm, developed
           by Bernard Rauscher.
        - `median`: Background noise is characterized by a median
          along the detector slow axis. Implementation is based on the
          `image1overf` algorithm, developed by Chris Willott.
    """

    class_alias = "clean_noise"

    spec = """
        algorithm = option('fft', 'median', default='fft')  # Cleaning algorithm
        mask_spectral_regions = boolean(default=True)  # Mask WCS-defined regions for spectral data
        n_integ = int(default=None)  # Number of integrations to use to create the mask, if not all
        n_sigma = float(default=5.0)  # Clipping level for outliers
        save_mask = boolean(default=False)  # Save the created mask
        user_mask = string(default=None)  # Path to user-supplied mask
        skip = boolean(default=True)  # By default, skip the step
    """

    def process(self, input):
        """
        Fit and subtract 1/f background noise from a ramp data set.

        Parameters
        ----------
        input : `~jwst.datamodels.RampModel`
            Input datamodel to be corrected

        n_sigma : float, optional
            Sigma clipping threshold to be used in detecting outliers in the image

        n_integ : int, optional
            Number of integrations to use to create the mask. Set to None
            to create a mask for each integration.

        save_mask : bool, optional
            Save the computed mask image

        user_mask : None, str, or `~jwst.datamodels.ImageModel`
            Optional user-supplied mask image; path to file or opened datamodel

        mask_spectral_regions : bool, optional
            Mask regions of the image defined by WCS bounding boxes for slits/slices

        Returns
        -------
        output_model : `~jwst.datamodels.RampModel`
            The 1/f corrected datamodel
        """

        # Open the input data model
        with datamodels.open(input) as input_model:

            result = clean_noise.do_correction(
                input_model, self.algorithm, self.mask_spectral_regions,
                self.n_sigma, self.n_integ, self.save_mask, self.user_mask)
            output_model, mask_model = result

            # Save the mask, if requested
            if self.save_mask and mask_model is not None:
                mask_path = self.make_output_path(basepath=input_model.meta.filename, suffix='mask')
                self.log.info(f"Saving mask file {mask_path}")
                mask_model.save(mask_path)
                mask_model.close()

        return output_model
