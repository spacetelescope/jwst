"""Master Background Pipeline for applying Master Background to NIRSpec Slit-like data"""
from . import nirspec_utils
from ..barshadow import barshadow_step
from .. import datamodels
from ..flatfield import flat_field_step
from ..pathloss import pathloss_step
from ..photom import photom_step
from ..stpipe.step import preserve_step_pars
from ..stpipe import Pipeline

__all__ = ['MasterBackgroundNrsSlitsStep']

# Step parameters to generally ignore when copying from the parent steps.
GLOBAL_PARS_TO_IGNORE = ['output_ext', 'output_file', 'output_use_model', 'output_use_index',
                         'inverse', 'pre_hooks', 'post_hooks', 'save_results', 'suffix']


class MasterBackgroundNrsSlitsStep(Pipeline):
    """Apply master background processing to NIRSpec Slit-like data

    For MOS, and ignoring FS, the calibration process needs to occur
    twice: Once to calibrate background slits and create a master background.
    Then a second time to calibrate science using the master background.

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

    spec = """
        force_subtract = boolean(default=False)  # Force subtracting master background
        save_background = boolean(default=False) # Save computed master background
        user_background = string(default=None)   # Path to user-supplied master background
        inverse = boolean(default=False)    # Invert the operation
        output_use_model = boolean(default=True)
    """

    # Define aliases to steps
    step_defs = {
        'flat_field': flat_field_step.FlatFieldStep,
        'pathloss': pathloss_step.PathLossStep,
        'barshadow': barshadow_step.BarShadowStep,
        'photom': photom_step.PhotomStep,
    }

    # No need to prefetch. This will have been done by the parent step.
    prefetch_references = False

    def process(self, data):
        """Compute and subtract a master background spectrum

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`
            The data to operate on.

        Attributes
        ----------
        correction_pars : dict
            The master background information from a previous invocation of the step.
            Keys are:

            - "masterbkg_1d": `~jwst.datamodels.CombinedSpecModel`
                The 1D version of the master background.

            - "masterbkg_2d": `~jwst.datamodels.MultiSlitModel`
                The 2D slit-based version of the master background.

        force_subtract : bool, optional
            Optional user-supplied flag that overrides step logic to force subtraction of the
            master background.
            Default is False, in which case the step logic determines if the calspec2 background step
            has already been applied and, if so, the master background step is skipped.
            If set to True, the step logic is bypassed and the master background is subtracted.

        save_background : bool, optional
            Save computed master background.

        user_background : None, string, or `~jwst.datamodels.CombinedSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel

        Returns
        -------
        result : `~jwst.datamodels.MultiSlitModel`
        """
        with datamodels.open(data) as data_model:
            # If some type of background processing had already been done. Abort.
            # UNLESS forcing is enacted.
            if not self.force_subtract and \
               'COMPLETE' in [data_model.meta.cal_step.back_sub, data_model.meta.cal_step.master_background]:
                self.log.info('Background subtraction has already occurred. Skipping.')
                self.record_step_status(data, 'master_background', False)
                return data

            if self.user_background:
                self.log.info(f'Calculating master background from user-supplied background {self.user_background}')
                user_background = datamodels.open(self.user_background)
                master_background, mb_multislit = self._calc_master_background(data_model, user_background)
            elif self.use_correction_pars:
                self.log.info('Using pre-calculated correction parameters.')
                master_background = self.correction_pars['masterbkg_1d']
                mb_multislit = self.correction_pars['masterbkg_2d']
            else:
                num_bkg, num_src = self._classify_slits(data_model)
                if num_bkg == 0:
                    self.log.warning('No background slits available for creating master background. Skipping')
                    self.record_step_status(data, 'master_background', False)
                    return data
                elif num_src == 0:
                    self.log.warning('No source slits for applying master background. Skipping')
                    self.record_step_status(data, 'master_background', False)
                    return data

                self.log.info('Calculating master background')
                master_background, mb_multislit = self._calc_master_background(data_model)

            # Check that a master background was actually determined.
            if master_background is None:
                self.log.info('No master background could be calculated. Skipping.')
                self.record_step_status(data, 'master_background', False)
                return data

            # Now apply the de-calibrated background to the original science
            result = nirspec_utils.apply_master_background(data_model, mb_multislit, inverse=self.inverse)

            # Mark as completed and setup return data
            self.record_step_status(result, 'master_background', True)
            self.correction_pars = {
                'masterbkg_1d': master_background,
                'masterbkg_2d': mb_multislit
            }
            if self.save_background:
                self.save_model(master_background, suffix='masterbg1d', force=True)
                self.save_model(mb_multislit, suffix='masterbg2d', force=True)

        return result

    def set_pars_from_parent(self):
        """Set substep parameters from the parents substeps"""
        if not self.parent:
            return

        steps = ['barshadow', 'flat_field', 'pathloss', 'photom']
        pars_to_ignore = {
            'barshadow': ['source_type'],
            'flat_field': ['save_interpolated_flat'],
            'pathloss': ['source_type'],
            'photom': ['source_type']
        }

        for step in steps:
            pars = getattr(self.parent, step).get_pars()
            for par in pars_to_ignore[step] + GLOBAL_PARS_TO_IGNORE:
                del pars[par]
            getattr(self, step).update_pars(pars)

    def _classify_slits(self, data):
        """Determine how many Slits are background and source types

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
            if "background" in slit.source_name:
                num_bkg+=1
            else:
                num_src+=1

        return num_bkg, num_src

    def _calc_master_background(self, data, user_background=None):
        """Calculate master background from background slits

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`
            The data to operate on.

        user_background : None, string, or `~jwst.datamodels.CombinedSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel

        Returns
        -------
        masterbkg_1d, masterbkg_2d : `~jwst.datamodels.CombinedSpecModel`, `~jwst.datamodels.MultiSlitModel`
            The master background in 1d and 2d, multislit formats.
            None is returned when a master background could not be determined.
        """

        # Since the parameters for the substeps are modified during processing,
        # wrap the processing in a context manager that restores all parameters.
        with preserve_step_pars(self):

            # When this step is called from another step/pipeliine,
            # retrieve the matching substep parameters from the parent.
            # This permits the substeps to perform similarly to what is
            # specified in the parent's substeps, such as skipping.
            # Any parameters that need be changed below are ignored.
            self.set_pars_from_parent()

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
            if user_background:
                self.log.debug(f'User background provided {user_background}')
                master_background = user_background
            else:
                self.log.debug('Calculating 1D master background')
                master_background = nirspec_utils.create_background_from_multislit(pre_calibrated)
            if master_background is None:
                self.log.debug('No master background could be calculated. Returning None')
                return None, None

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

        return master_background, mb_multislit
