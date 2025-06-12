from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise import clean_flicker_noise
from jwst.stpipe import Step

__all__ = ["NSCleanStep"]


class NSCleanStep(Step):
    """
    Perform 1/f noise correction.

    NSCleanStep: This step performs 1/f noise correction ("cleaning")
    of NIRSpec images, using the "NSClean" method.

    NOTE: This step is a deprecated alias to the ``clean_flicker_noise`` step.

    Attributes
    ----------
    fit_method : str, optional
        The background fit algorithm to use.  Options are 'fft' and 'median';
        'fft' performs the original NSClean implementation.
    fit_by_channel : bool, optional
        If set, flicker noise is fit independently for each detector channel.
        Ignored for subarray data and for `fit_method` = 'fft'.
    background_method : {'median', 'model', None}
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is fit with a low-resolution model via
        `~photutils.background.Background2D`.  If None, the background
        value is 0.0.
    background_box_size : tuple of int, optional
        Box size for the data grid used by `Background2D` when
        `background_method` = 'model'. For best results, use a box size
        that evenly divides the input image shape.
    mask_spectral_regions : bool, optional
        Mask regions of the image defined by WCS bounding boxes for slits/slices.
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
        self.log.warning(message)

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
            output_model.meta.cal_step.nsclean = status

        return output_model
