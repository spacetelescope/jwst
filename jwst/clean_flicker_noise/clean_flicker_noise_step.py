from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import clean_flicker_noise

__all__ = ["CleanFlickerNoiseStep"]


class CleanFlickerNoiseStep(Step):
    """
    Perform flicker noise correction.

    Input data is expected to be a ramp file (RampModel), in between
    jump and ramp fitting steps, or a rate file (ImageModel or CubeModel).

    Correction algorithms implemented are:
        - `fft`: Background noise is fit in frequency space.
           Implementation is based on the NSClean algorithm, developed
           by Bernard Rauscher.
        - `median`: Background noise is characterized by a median
          along the detector slow axis. Implementation is based on the
          `image1overf` algorithm, developed by Chris Willott.

    Attributes
    ----------
    fit_method : str, optional
        The noise fitting algorithm to use.  Options are 'fft' and 'median'.
    fit_by_channel : bool, optional
        If set, flicker noise is fit independently for each detector channel.
        Ignored for MIRI, for subarray data, and for `fit_method` = 'fft'
    background_method : {'median', 'model', None}
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is fit with a low-resolution model via
        `~photutils.background.Background2D`.  If None, the background
        value is 0.0.
    background_box_size : tuple of int or None, optional
        Box size for the data grid used by `Background2D` when
        `background_method` = 'model'. For best results, use a box size
        that evenly divides the input image shape.  If None, the largest
        value between 1 and 32 that evenly divides the image dimension
        is used.
    mask_science_regions : bool, optional
        For NIRSpec, mask regions of the image defined by WCS bounding
        boxes for slits/slices, as well as any regions known to be
        affected by failed-open MSA shutters.  For MIRI imaging, mask
        regions of the detector not used for science.
    apply_flat_field : bool, optional
        If set, images are flat-corrected prior to fitting background
        and noise levels.  A full-frame flat field image
        (reference type FLAT) is required. For modes that do not provide
        FLAT files via CRDS, including all NIRSpec spectral modes, a manually
        generated override flat is required to enable this option.
        Use the `override_flat` parameter to provide an alternate flat image
        as needed.
    n_sigma : float, optional
        Sigma clipping threshold to be used in detecting outliers in the image.
    fit_histogram : bool, optional
        If set, the 'sigma' used with `n_sigma` for clipping outliers
        is derived from a Gaussian fit to a histogram of values.
        Otherwise, a simple iterative sigma clipping is performed.
    single_mask : bool, optional
        If set, a single mask will be created, regardless of
        the number of input integrations. Otherwise, the mask will
        be a 3D cube, with one plane for each integration.
    user_mask : None, str, or `~jwst.datamodels.ImageModel`
        Optional user-supplied mask image; path to file or opened datamodel.
    save_mask : bool, optional
        Save the computed mask image.
    save_background : bool, optional
        Save the computed background image.
    save_noise : bool, optional
        Save the computed noise image.
    """

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
                    self.log.warning(
                        f"Flat correction is not available for "
                        f"exposure type {exp_type} without a user-"
                        f"supplied flat."
                    )
                    flat_filename = None
                else:
                    self.log.info(f"Using FLAT reference file: {flat_filename}")

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
                self.log.info(f"Saving mask file {mask_path}")
                mask_model.save(mask_path)
                mask_model.close()

            # Save the background, if requested
            if self.save_background and background_model is not None:
                bg_path = self.make_output_path(
                    basepath=input_model.meta.filename, suffix="flicker_bkg"
                )
                self.log.info(f"Saving background file {bg_path}")
                background_model.save(bg_path)
                background_model.close()

            # Save the noise, if requested
            if self.save_noise and noise_model is not None:
                noise_path = self.make_output_path(
                    basepath=input_model.meta.filename, suffix="flicker_noise"
                )
                self.log.info(f"Saving noise file {noise_path}")
                noise_model.save(noise_path)
                noise_model.close()

            # Set the step completion status
            output_model.meta.cal_step.clean_flicker_noise = status

        return output_model
