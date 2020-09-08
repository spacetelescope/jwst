"""Master Background Pipeline for applying Master Background to NIRSpec Slit-like data"""

from ..barshadow import barshadow_step
from ..flatfield import flat_field_step
from ..master_background import nirspec_utils
from ..pathloss import pathloss_step
from ..photom import photom_step
from ..stpipe import Pipeline


class MasterBackgroundNRSSlitsPipe(Pipeline):
    """Apply master background processing to NIRSpec Slit-like data

    For MOS, and ignoring FS, the calibration process needs to occur
    twice: Once to calibrate background slits and create a master background.
    Then a second time to calibrate science using the master background.

    Notes
    -----
    The algorithm is as follows:

    - Calibrate all slits
      - For each step:
        - Force the source type to be extended source for all slits.
        - Return the correction array used.
    - Create the 1D master background
    - For each slit
      - Expand out the 1D master background to match the 2D wavelength grid of the slit
      - Reverse-calibrate the 2D background, using the correction arrays calculated above.
      - Subtract the background from the input slit data
    """

    spec = """
        user_background = string(default=None)   # Path to user-supplied master background
        save_background = boolean(default=False) # Save computed master background
        force_subtract = boolean(default=False)  # Force subtracting master background
        output_use_model = boolean(default=True)
    """

    # Define aliases to steps
    step_defs = {
        'flat_field': flat_field_step.FlatFieldStep,
        'pathloss': pathloss_step.PathLossStep,
        'barshadow': barshadow_step.BarShadowStep,
        'photom': photom_step.PhotomStep,
    }

    def process(self, data):
        """Compute and subtract a master background spectrum

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`

        user_background : None, string, or `~jwst.datamodels.MultiSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel

        save_background : bool, optional
            Save computed master background.

        force_subtract : bool, optional
            Optional user-supplied flag that overrides step logic to force subtraction of the
            master background.
            Default is False, in which case the step logic determines if the calspec2 background step
            has already been applied and, if so, the master background step is skipped.
            If set to True, the step logic is bypassed and the master background is subtracted.

        Returns
        -------
        result : `~jwst.datamodels.MultiSlitModel`
        """

        # If some type of background processing had already been done. Abort.
        # UNLESS forcing is enacted.
        if not self.force_subtract and \
           'COMPLETE' in [data.meta.cal_step.back_sub, data.meta.cal_step.master_background]:
            self.log.info('Background subtraction has already occurred. Skipping')
            data.meta.cal_step.master_background = 'SKIP'
            return data

        # First pass: just do the calibration to determine the correction
        # arrays. However, force all slits to be processed as extended sources.
        self.pathloss.source_type = 'EXTENDED'
        self.barshadow.source_type = 'EXTENDED'
        self.photom.source_type = 'EXTENDED'

        pre_calibrated = self.flat_field(data)
        pre_calibrated = self.pathloss(pre_calibrated)
        pre_calibrated = self.barshadow(pre_calibrated)
        pre_calibrated = self.photom(pre_calibrated)

        # Create the 1D, fully calibrated master background.
        master_background = nirspec_utils.create_background_from_multislit(pre_calibrated)
        if master_background is None:
            return data

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

        mb_multislit = self.photom(mb_multislit)
        mb_multislit = self.barshadow(mb_multislit)
        mb_multislit = self.pathloss(mb_multislit)
        mb_multislit = self.flat_field(mb_multislit)

        # Now apply the de-calibrated background to the original science
        # At this point, should just be a slit-to-slit subtraction operation.
        result = nirspec_utils.apply_master_background(data, mb_multislit)

        # Mark as completed.
        result.meta.cal_step.master_background = 'COMPLETE'

        return result
