import logging

from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise import clean_flicker_noise
from jwst.stpipe import Step

__all__ = ["CleanFlickerNoiseStep"]

log = logging.getLogger(__name__)


class CleanFlickerNoiseStep(Step):
    """Perform flicker noise correction."""

    class_alias = "clean_flicker_noise"

    spec = """
        fit_method = option('fft', 'median', default='median')  # Noise fitting algorithm
        fit_by_channel = boolean(default=False)  # Fit noise separately by amplifier (NIR only)
        background_method = option('median', 'model', None, default='median') # Background fit
        background_box_size = int_list(min=2, max=2, default=None)  # Background box size
        mask_science_regions = boolean(default=False)  # Mask known science regions
        apply_flat_field = boolean(default=False)  # Apply a flat correction before fitting
        n_sigma = float(default=2.0)  # Clipping level for non-background signal
        fit_histogram = boolean(default=False)  # Fit a value histogram to derive sigma
        single_mask = boolean(default=True)  # Make a single mask for all integrations
        user_mask = string(default=None)  # Path to user-supplied mask
        save_mask = boolean(default=False)  # Save the created mask
        save_background = boolean(default=False)  # Save the fit background
        save_noise = boolean(default=False)  # Save the fit noise
        skip = boolean(default=True)  # By default, skip the step
    """  # noqa: E501

    reference_file_types = ["flat"]

    def process(self, input_data):
        """
        Fit and subtract 1/f background noise from a ramp data set.

        Input data is expected to be a ramp file (RampModel), in between
        jump and ramp fitting steps, or a rate file (ImageModel or CubeModel).

        Correction algorithms implemented are:

            - "fft": Background noise is fit in frequency space.
               Implementation is based on the NSClean algorithm, developed
               by Bernard Rauscher.
            - "median": Background noise is characterized by a median
              along the detector slow axis. Implementation is based on the
              "image1overf" algorithm, developed by Chris Willott.

        Parameters
        ----------
        input_data : DataModel
            Input datamodel to be corrected

        Returns
        -------
        output_model : DataModel
            The flicker noise corrected datamodel
        """
        # Open the input data model
        with datamodels.open(input_data) as input_model:
            flat_filename = None
            if self.apply_flat_field:
                flat_filename = self.get_reference_file(input_model, "flat")
                exp_type = input_model.meta.exposure.type
                if flat_filename == "N/A":
                    log.warning(
                        f"Flat correction is not available for "
                        f"exposure type {exp_type} without a user-"
                        f"supplied flat."
                    )
                    flat_filename = None
                else:
                    log.info(f"Using FLAT reference file: {flat_filename}")

            result = clean_flicker_noise.do_correction(
                input_model,
                input_dir=self.input_dir,
                fit_method=self.fit_method,
                fit_by_channel=self.fit_by_channel,
                background_method=self.background_method,
                background_box_size=self.background_box_size,
                mask_science_regions=self.mask_science_regions,
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
            output_model.meta.cal_step.clean_flicker_noise = status

        return output_model
