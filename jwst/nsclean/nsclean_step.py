import logging

from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise import clean_flicker_noise
from jwst.stpipe import Step

__all__ = ["NSCleanStep"]

log = logging.getLogger(__name__)


class NSCleanStep(Step):
    """
    Perform 1/f noise correction.

    NOTE: This step is a deprecated alias to the ``clean_flicker_noise`` step.
    """

    class_alias = "nsclean"

    spec = """
        fit_method = option('fft', 'median', default='fft')  # Noise fitting algorithm
        fit_by_channel = boolean(default=False)  # Fit noise separately by amplifier
        background_method = option('median', 'model', None, default=None) # Background fit
        background_box_size = int_list(min=2, max=2, default=None)  # Background box size
        mask_spectral_regions = boolean(default=True)  # Mask WCS-defined spectral regions
        n_sigma = float(default=5.0)  # Clipping level for outliers
        fit_histogram = boolean(default=False)  # Fit a value histogram to derive sigma
        single_mask = boolean(default=False)  # Make a single mask for all integrations
        user_mask = string(default=None)  # Path to user-supplied mask
        save_mask = boolean(default=False)  # Save the created mask
        save_background = boolean(default=False)  # Save the fit background
        save_noise = boolean(default=False)  # Save the fit noise
        skip = boolean(default=True)  # By default, skip the step
    """  # noqa: E501

    def process(self, input_data):
        """
        Fit and subtract 1/f background noise from a NIRSpec image.

        Parameters
        ----------
        input_data : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`
            Input datamodel to be corrected.

        Returns
        -------
        output_model : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`
            The 1/f corrected datamodel.
        """
        message = (
            "The 'nsclean' step is a deprecated alias to 'clean_flicker_noise' "
            "and will be removed in future builds."
        )
        log.warning(message)

        # Open the input data model
        with datamodels.open(input_data) as input_model:
            # clean_flicker_noise allows flat handling, but NIRSpec is
            # not supported via CRDS files, since it does not have a full
            # frame flat.  Since this step is for NIRSpec only, don't allow
            # flat handling here.
            flat_filename = None

            # Do the NSClean correction
            result = clean_flicker_noise.do_correction(
                input_model,
                input_dir=self.input_dir,
                fit_method=self.fit_method,
                fit_by_channel=self.fit_by_channel,
                background_method=self.background_method,
                background_box_size=self.background_box_size,
                mask_science_regions=self.mask_spectral_regions,
                flat_filename=flat_filename,
                n_sigma=self.n_sigma,
                fit_histogram=self.fit_histogram,
                single_mask=self.single_mask,
                user_mask=self.user_mask,
                save_mask=self.save_mask,
                save_background=self.save_background,
                save_noise=self.save_noise,
            )
            output_model, mask_model, background_model, noise_model, status = result

            # Save the mask, if requested
            if self.save_mask and mask_model is not None:
                mask_path = self.make_output_path(basepath=input_model.meta.filename, suffix="mask")
                log.info(f"Saving mask file {mask_path}")
                mask_model.save(mask_path)
                mask_model.close()

            # Save the background, if requested
            if self.save_background and background_model is not None:
                bg_path = self.make_output_path(
                    basepath=input_model.meta.filename, suffix="flicker_bkg"
                )
                log.info(f"Saving background file {bg_path}")
                background_model.save(bg_path)
                background_model.close()

            # Save the noise, if requested
            if self.save_noise and noise_model is not None:
                noise_path = self.make_output_path(
                    basepath=input_model.meta.filename, suffix="flicker_noise"
                )
                log.info(f"Saving noise file {noise_path}")
                noise_model.save(noise_path)
                noise_model.close()

            # Set the step completion status
            output_model.meta.cal_step.nsclean = status

        return output_model
