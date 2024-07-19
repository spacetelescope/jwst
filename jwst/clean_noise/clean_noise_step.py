from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import clean_noise

__all__ = ["CleanNoiseStep"]


class CleanNoiseStep(Step):
    """
    CleanNoiseStep: This step performs 1/f noise correction ("cleaning").

    Input data is expected to be a ramp file (RampModel), in between
    jump and ramp fitting steps, or a rate file (ImageModel or CubeModel).

    Correction algorithms implemented are:
        - `fft`: Background noise is fit in frequency space.
           Implementation is based on the NSClean algorithm, developed
           by Bernard Rauscher.
        - `median`: Background noise is characterized by a median
          along the detector slow axis. Implementation is based on the
          `image1overf` algorithm, developed by Chris Willott.
    """

    class_alias = "clean_noise"

    spec = """
        fit_method = option('fft', 'median', default='median')  # Noise fitting algorithm
        background_method = option('median', 'model', None, default='median')
        mask_spectral_regions = boolean(default=False)  # Mask WCS-defined regions for spectral data
        single_mask = boolean(default=False)  # Make a single mask for all integrations
        n_sigma = float(default=5.0)  # Clipping level for outliers
        save_mask = boolean(default=False)  # Save the created mask
        user_mask = string(default=None)  # Path to user-supplied mask
        use_diff = boolean(default=False)  # Correct group diffs instead of group images
        skip = boolean(default=True)  # By default, skip the step        
    """

    def process(self, input):
        """
        Fit and subtract 1/f background noise from a ramp data set.

        Parameters
        ----------
        input : DataModel
            Input datamodel to be corrected

        fit_method : str, optional
            The background fit algorithm to use.  Options are 'fft' and 'median'.

        background_method : {'median', 'model', None}, optional
            If 'median', the preliminary background to remove and restore
            is a simple median of the background data.  If 'model', the
            background data is modeled with a median filter with a 5x5
            pixel kernel.  If None, the background value is 0.0.

        mask_spectral_regions : bool, optional
            Mask regions of the image defined by WCS bounding boxes for slits/slices.
            Ignored for imaging data.

        n_sigma : float, optional
            Sigma clipping threshold to be used in detecting outliers in the image.

        single_mask : bool, optional
            If set, a single mask will be created, regardless of
            the number of input integrations. Otherwise, the mask will
            be a 3D cube, with one plane for each integration.

        save_mask : bool, optional
            Save the computed mask image.

        user_mask : None, str, or DataModel
            Optional user-supplied mask image; path to file or opened datamodel.
            May be a 2D image (ImageModel) or 3D cube (CubeModel), matching the
            number of integrations in the data.

        use_diff : bool, optional
            If set, and the input is ramp data, correction is performed
            on diffs between group images.  Otherwise, correction is
            performed directly on the group image.

        Returns
        -------
        output_model : DataModel
            The 1/f corrected datamodel
        """

        # Open the input data model
        with datamodels.open(input) as input_model:

            result = clean_noise.do_correction(
                input_model, self.fit_method, self.background_method,
                self.mask_spectral_regions,
                self.n_sigma, self.single_mask, self.save_mask, self.user_mask,
                self.use_diff)
            output_model, mask_model, status = result

            # Save the mask, if requested
            if self.save_mask and mask_model is not None:
                mask_path = self.make_output_path(basepath=input_model.meta.filename, suffix='mask')
                self.log.info(f"Saving mask file {mask_path}")
                mask_model.save(mask_path)
                mask_model.close()

            # Set the step completion status
            output_model.meta.cal_step.clean_noise = status

        return output_model
