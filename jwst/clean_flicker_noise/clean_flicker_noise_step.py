import logging

from jwst.clean_flicker_noise import autoparam, clean_flicker_noise
from jwst.stpipe import Step

__all__ = ["CleanFlickerNoiseStep"]

log = logging.getLogger(__name__)


class CleanFlickerNoiseStep(Step):
    """Perform flicker noise correction."""

    class_alias = "clean_flicker_noise"

    spec = """
        autoparam = boolean(default=False) # Automatically select some fit and background parameters for the input data.
        fit_method = option('fft', 'median', default='median')  # Noise fitting algorithm.
        fit_by_channel = boolean(default=False)  # Fit noise separately by amplifier (NIR only).
        background_method = option('median', 'model', None, default='median') # Background fit.
        background_box_size = int_list(min=2, max=2, default=None)  # Background box size.
        mask_science_regions = boolean(default=False)  # Mask known science regions.
        apply_flat_field = boolean(default=False)  # Apply a flat correction before fitting.
        n_sigma = float(default=2.0)  # Clipping level for non-background signal.
        fit_histogram = boolean(default=False)  # Fit a value histogram to derive sigma.
        single_mask = boolean(default=True)  # Make a single mask for all integrations.
        user_mask = string(default=None)  # Path to user-supplied mask
        save_mask = boolean(default=False)  # Save the created mask
        save_background = boolean(default=False)  # Save the fit background
        save_noise = boolean(default=False)  # Save the fit noise
        skip = boolean(default=True)  # By default, skip the step.
    """  # noqa: E501

    reference_file_types = ["flat"]

    def _set_auto_parameters(self, input_model):
        """
        Override fit parameters from input data characteristics if possible.

        For any supported exposure type, the data is inspected and fit
        parameters are determined accordingly.  Any parameters specified
        by the `autoparam` algorithm are directly overridden, ignoring
        user input.

        Parameters
        ----------
        input_model : DataModel
            Input datamodel to be corrected.
        """
        exp_type = input_model.meta.exposure.type
        found_exptype = True
        if exp_type == "NIS_IMAGE":
            flat_filename = self.get_reference_file(input_model, "flat")
            override_parameters = autoparam.niriss_image_parameters(input_model, flat_filename)
        elif exp_type == "NRC_IMAGE":
            flat_filename = self.get_reference_file(input_model, "flat")
            override_parameters = autoparam.nircam_image_parameters(input_model, flat_filename)
        else:
            override_parameters = None
            found_exptype = False

        if override_parameters is not None:
            log.info(f"Auto parameters set for {exp_type}:")
            for param, value in override_parameters.items():
                log.info(f"  {param}: {value}")
                setattr(self, param, value)
        else:
            if found_exptype:
                log.warning("Auto parameter setting failed.")
            else:
                log.warning(f"Auto parameters are not available for exposure type {exp_type}")
            log.info("Using input parameters as provided; no overrides applied.")

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
        input_data : str or `~stdatamodels.jwst.datamodels.RampModel` \
                     or `~stdatamodels.jwst.datamodels.ImageModel` \
                     or `~stdatamodels.jwst.datamodels.CubeModel`
            Filename or input datamodel to be corrected.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.RampModel` \
                       or `~stdatamodels.jwst.datamodels.ImageModel` \
                       or `~stdatamodels.jwst.datamodels.CubeModel`
            The flicker noise corrected datamodel, matching the type of
            the input.
        """
        # Open the input data model
        # Note: this step can be run either in stage 1 on ramp-type
        # data or in stage 2 on rate-type data, so it's important
        # not to provide a specific datamodel to prepare_output.
        output_model = self.prepare_output(input_data)

        # Assign fit and background parameters appropriate
        # to the input data if desired
        if self.autoparam:
            self._set_auto_parameters(output_model)

        flat_filename = None
        if self.apply_flat_field:
            flat_filename = self.get_reference_file(output_model, "flat")
            exp_type = output_model.meta.exposure.type
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
            output_model,
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
            mask_path = self.make_output_path(basepath=output_model.meta.filename, suffix="mask")
            log.info(f"Saving mask file {mask_path}")
            mask_model.save(mask_path)
            mask_model.close()

        # Save the background, if requested
        if self.save_background and background_model is not None:
            bg_path = self.make_output_path(
                basepath=output_model.meta.filename, suffix="flicker_bkg"
            )
            log.info(f"Saving background file {bg_path}")
            background_model.save(bg_path)
            background_model.close()

        # Save the noise, if requested
        if self.save_noise and noise_model is not None:
            noise_path = self.make_output_path(
                basepath=output_model.meta.filename, suffix="flicker_noise"
            )
            log.info(f"Saving noise file {noise_path}")
            noise_model.save(noise_path)
            noise_model.close()

        # Set the step completion status
        output_model.meta.cal_step.clean_flicker_noise = status

        return output_model
