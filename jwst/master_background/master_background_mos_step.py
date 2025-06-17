"""Apply master background corrections to NIRSpec MOS data."""

from stpipe.step import preserve_step_pars
from jwst.stpipe import record_step_status

from stdatamodels.jwst import datamodels

from . import nirspec_utils
from jwst.barshadow import barshadow_step
from jwst.flatfield import flat_field_step
from jwst.pathloss import pathloss_step
from jwst.photom import photom_step
from jwst.pixel_replace import pixel_replace_step
from jwst.resample import resample_spec_step
from jwst.extract_1d import extract_1d_step
from jwst.stpipe import Pipeline

__all__ = ["MasterBackgroundMosStep"]

# Step parameters to generally ignore when copying from the parent steps.
GLOBAL_PARS_TO_IGNORE = [
    "output_ext",
    "output_file",
    "output_use_model",
    "output_use_index",
    "inverse",
    "pre_hooks",
    "post_hooks",
    "save_results",
    "suffix",
]


class MasterBackgroundMosStep(Pipeline):
    """
    Apply master background processing to NIRSpec MOS data.

    For MOS, and ignoring FS, the calibration process needs to occur
    twice: Once to calibrate background slits and create a master background.
    Then a second time to calibrate science using the master background.

    Attributes
    ----------
    correction_pars : dict
        The master background information from a previous invocation of the step.
        Keys are:

        - "masterbkg_1d": `~jwst.datamodels.CombinedSpecModel`
            The 1D version of the master background.
        - "masterbkg_2d": `~jwst.datamodels.MultiSlitModel`
            The 2D slit-based version of the master background.

    sigma_clip : None or float
        Optional factor for sigma clipping outliers when combining background spectra.
    median_kernel : int
        Optional user-supplied kernel with which to moving-median boxcar filter
        the master background spectrum.  Must be an odd integer; even integers will be
        rounded down to the nearest odd integer.
    force_subtract : bool
        Optional user-supplied flag that overrides step logic to force subtraction of the
        master background.
        Default is False, in which case the step logic determines if the calspec2 background step
        has already been applied and, if so, the master background step is skipped.
        If set to True, the step logic is bypassed and the master background is subtracted.
    save_background : bool
        Save computed master background.
    user_background : None, str, or `~jwst.datamodels.CombinedSpecModel`
        Optional user-supplied master background 1D spectrum, path to file
        or opened datamodel

    Notes
    -----
    The algorithm is as follows

    - Calibrate all slits

      - For each step

        - Force the source type to be extended source for all slits.
        - Return the correction array used.

    - Create the 1D master background
    - For each slit

      - Expand out the 1D master background to match the 2D wavelength grid of the slit
      - Reverse-calibrate the 2D background, using the correction arrays calculated above.
      - Subtract the background from the input slit data
    """

    class_alias = "master_background_mos"

    spec = """
        sigma_clip = float(default=3) # Factor for clipping outliers when combining spectra
        median_kernel = integer(default=1) # Moving-median boxcar size to filter background spec
        force_subtract = boolean(default=False)  # Force subtracting master background
        save_background = boolean(default=False) # Save computed master background
        user_background = string(default=None)   # Path to user-supplied master background
        inverse = boolean(default=False)    # Invert the operation
        output_use_model = boolean(default=True)
    """  # noqa: E501

    # Define aliases to steps
    step_defs = {
        "flat_field": flat_field_step.FlatFieldStep,
        "pathloss": pathloss_step.PathLossStep,
        "barshadow": barshadow_step.BarShadowStep,
        "photom": photom_step.PhotomStep,
        "pixel_replace": pixel_replace_step.PixelReplaceStep,
        "resample_spec": resample_spec_step.ResampleSpecStep,
        "extract_1d": extract_1d_step.Extract1dStep,
    }

    # No need to prefetch. This will have been done by the parent step.
    prefetch_references = False

    def process(self, data):
        """
        Compute and subtract a master background spectrum.

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`
            The data to operate on.

        Returns
        -------
        result : `~jwst.datamodels.MultiSlitModel`
            The background corrected data.
        """
        with datamodels.open(data) as data_model:
            # If some type of background processing had already been done. Abort.
            # UNLESS forcing is enacted.
            if not self.force_subtract and "COMPLETE" in [
                data_model.meta.cal_step.bkg_subtract,
                data_model.meta.cal_step.master_background,
            ]:
                self.log.info("Background subtraction has already occurred. Skipping.")
                record_step_status(data, "master_background", success=False)
                return data

            if self.user_background:
                self.log.info(
                    "Calculating master background from "
                    f"user-supplied background {self.user_background}"
                )
                user_background = datamodels.open(self.user_background)
                master_background, mb_multislit, bkg_x1d_spectra = self._calc_master_background(
                    data_model, user_background
                )
            elif self.use_correction_pars:
                self.log.info("Using pre-calculated correction parameters.")
                master_background = self.correction_pars["masterbkg_1d"]
                mb_multislit = self.correction_pars["masterbkg_2d"]
            else:
                num_bkg, num_src = self._classify_slits(data_model)
                if num_bkg == 0:
                    self.log.warning(
                        "No background slits available for creating master background. Skipping"
                    )
                    record_step_status(data, "master_background", False)
                    return data
                elif num_src == 0:
                    self.log.warning("No source slits for applying master background. Skipping")
                    record_step_status(data, "master_background", False)
                    return data

                self.log.info("Calculating master background")
                master_background, mb_multislit, bkg_x1d_spectra = self._calc_master_background(
                    data_model, sigma_clip=self.sigma_clip, median_kernel=self.median_kernel
                )

            # Check that a master background was actually determined.
            if master_background is None:
                self.log.info("No master background could be calculated. Skipping.")
                record_step_status(data, "master_background", False)
                return data

            # Now apply the de-calibrated background to the original science
            result = nirspec_utils.apply_master_background(
                data_model, mb_multislit, inverse=self.inverse
            )

            # Mark as completed and setup return data
            record_step_status(result, "master_background", True)
            self.correction_pars = {"masterbkg_1d": master_background, "masterbkg_2d": mb_multislit}
            if self.save_background:
                self.save_model(master_background, suffix="masterbg1d", force=True)
                self.save_model(mb_multislit, suffix="masterbg2d", force=True)
                if bkg_x1d_spectra is not None:
                    self.save_model(bkg_x1d_spectra, suffix="bkgx1d", force=True)

        return result

    def set_pars_from_parent(self):
        """Set substep parameters from the parents substeps when needed."""
        if not self.parent:
            return

        steps = ["barshadow", "flat_field", "pathloss", "photom"]
        pars_to_ignore = {
            "barshadow": ["source_type"],
            "flat_field": ["save_interpolated_flat"],
            "pathloss": ["source_type"],
            "photom": ["source_type"],
        }

        for step in steps:
            pars = getattr(self.parent, step).get_pars()
            for par in pars_to_ignore[step] + GLOBAL_PARS_TO_IGNORE:
                del pars[par]
            getattr(self, step).update_pars(pars)

    def _extend_bg_slits(self, pre_calibrated):
        # Copy dedicated background slitlets to a temporary model
        bkg_model = datamodels.MultiSlitModel()
        bkg_model.update(pre_calibrated)
        slits = []
        for slit in pre_calibrated.slits:
            if nirspec_utils.is_background_msa_slit(slit):
                self.log.info(f"Using background slitlet {slit.source_name}")
                slits.append(slit)
        if len(slits) == 0:
            self.log.warning("No background slitlets found; skipping master bkg correction")
            return None
        bkg_model.slits.extend(slits)
        return bkg_model

    def _classify_slits(self, data):
        """
        Determine how many slits are background and source types.

        Parameters
        ----------
        data : ~jwst.datamodels.MultiSlitModel`
            The data to operate on.

        Returns
        -------
        num_bkg, num_src : int, int
            The number of background slits and the number of source slits.
        """
        # Loop over all the Slit instances in the input data model and
        # count how many are background vs source.
        num_bkg = num_src = 0
        for slit in data.slits:
            if nirspec_utils.is_background_msa_slit(slit):
                num_bkg += 1
            else:
                num_src += 1

        return num_bkg, num_src

    def _calc_master_background(self, data, user_background=None, sigma_clip=3, median_kernel=1):
        """
        Calculate master background from background slits.

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`
            The data to operate on.
        user_background : None, str, or `~jwst.datamodels.CombinedSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel
        sigma_clip : None or float
            Optional factor for sigma clipping outliers when combining background spectra.
        median_kernel : int, optional
            Optional user-supplied kernel with which to moving-median boxcar
            filter the master background spectrum.  Must be an odd integer; even
            integers will be rounded down to the nearest odd integer.

        Returns
        -------
        masterbkg_1d : `~jwst.datamodels.CombinedSpecModel`
            The master background in 1d multislit format.
            None is returned when a master background could not be determined.
        masterbkg_2d : `~jwst.datamodels.MultiSlitModel`
            The master background in 2d, multislit format.
            None is returned when a master background could not be determined.
        bkg_x1d_spectra : `~jwst.datamodels.MultiSlitModel`
            The 1D extracted background spectra used to determine the master background.
            Returns None when a user_background is provided.
        """
        # Since the parameters for the substeps are modified during processing,
        # wrap the processing in a context manager that restores all parameters.
        with preserve_step_pars(self):
            # When this step is called from another step/pipeline,
            # retrieve the matching substep parameters from the parent.
            # This permits the substeps to perform similarly to what is
            # specified in the parent's substeps, such as skipping.
            # Any parameters that need be changed below are ignored.
            self.set_pars_from_parent()

            # First pass: just do the calibration to determine the correction
            # arrays. However, force all slits to be processed as extended sources.
            self.pathloss.source_type = "EXTENDED"
            self.barshadow.source_type = "EXTENDED"
            self.photom.source_type = "EXTENDED"

            pre_calibrated = self.flat_field.run(data)
            pre_calibrated = self.pathloss.run(pre_calibrated)
            pre_calibrated = self.barshadow.run(pre_calibrated)
            pre_calibrated = self.photom.run(pre_calibrated)

            # Create the 1D, fully calibrated master background.
            if user_background:
                self.log.debug(f"User background provided {user_background}")
                master_background = user_background
                bkg_x1d_spectra = None
            else:
                self.log.info("Creating MOS master background from background slitlets")
                bkg_model = self._extend_bg_slits(pre_calibrated)
                if bkg_model is not None:
                    bkg_model = self.pixel_replace.run(bkg_model)
                    bkg_model = self.resample_spec.run(bkg_model)
                    bkg_x1d_spectra = self.extract_1d.run(bkg_model)
                    master_background = nirspec_utils.create_background_from_multispec(
                        bkg_x1d_spectra, sigma_clip=sigma_clip, median_kernel=median_kernel
                    )
                else:
                    master_background = None
            if master_background is None:
                self.log.debug("No master background could be calculated. Returning None")
                return None, None, None

            # Now decalibrate the master background for each individual science slit.
            # First step is to map the master background into a MultiSlitModel
            # where the science slits are replaced by the master background.
            # Here the broadcasting from 1D to 2D need also occur.
            mb_multislit = nirspec_utils.map_to_science_slits(pre_calibrated, master_background)

            # Now that the master background is pretending to be science,
            # walk backwards through the steps to uncalibrate, using the
            # calibration factors carried from `pre_calibrated`.
            self.photom.use_correction_pars = True
            self.photom.inverse = True
            self.barshadow.use_correction_pars = True
            self.barshadow.inverse = True
            self.pathloss.use_correction_pars = True
            self.pathloss.inverse = True
            self.flat_field.use_correction_pars = True
            self.flat_field.inverse = True

            mb_multislit = self.photom.run(mb_multislit)
            mb_multislit = self.barshadow.run(mb_multislit)
            mb_multislit = self.pathloss.run(mb_multislit)
            mb_multislit = self.flat_field.run(mb_multislit)

        return master_background, mb_multislit, bkg_x1d_spectra
