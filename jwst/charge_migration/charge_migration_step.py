"""Detect and flag charge migration."""

import logging

from stdatamodels.jwst import datamodels

from jwst.charge_migration.charge_migration import correct_charge_migration, flag_charge_migration
from jwst.stpipe import Step

log = logging.getLogger(__name__)

__all__ = ["ChargeMigrationStep"]


class ChargeMigrationStep(Step):
    """Flag or correct groups exhibiting significant charge migration."""

    class_alias = "charge_migration"

    spec = """
        signal_threshold = float(default=-1)  # Threshold value for flagging excessive charge loss.
        correction_threshold = float(default=-1)  # Threshold value to correct for charge migration.
        mask_pad = integer(default=2)  # Number of pixels to fit on either side of a pixel with charge loss.
        snr_limit = float(default=3.0)  # Correction fit significance limit, as a factor of the read noise
        rms_limit = float(default=0.45)  # Correction fit quality limit, as a fraction of the initial RMS
        save_fit_info = boolean(default=False)  # Save a table of information about fit corrections
        skip = boolean(default=True)  # By default, skip the step.
    """  # noqa: E501

    reference_file_types = ["readnoise"]

    def process(self, step_input):
        """
        Flag or correct for charge migration effects.

        Parameters
        ----------
        step_input : `~stdatamodels.jwst.datamodels.RampModel`
            The ramp model on which to detect jumps.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            The flagged ramp model.
        """
        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        if result.data.shape[1] < 3:  # skip step if only 1 or 2 groups/integration
            log.info("Too few groups per integration; skipping charge_migration")
            result.meta.cal_step.charge_migration = "SKIPPED"
            return result

        if self.signal_threshold < 0 and self.correction_threshold < 0:
            log.info("All thresholds negative; skipping charge_migration")
            result.meta.cal_step.charge_migration = "SKIPPED"
            return result

        if self.signal_threshold >= 0:
            log.info("Flagging pixels with significant charge loss")
            result = flag_charge_migration(result, self.signal_threshold)

        if self.correction_threshold >= 0:
            log.info("Fitting and correcting for the brighter-fatter effect (BFE)")
            log.warning(
                "This process is currently a prototype: it is not yet recommended for general use."
            )
            readnoise_filename = self.get_reference_file(result, "readnoise")
            log.info("Using READNOISE reference file: %s", readnoise_filename)

            with datamodels.ReadnoiseModel(readnoise_filename) as rnoise:
                output = correct_charge_migration(
                    result,
                    rnoise,
                    correction_threshold=self.correction_threshold,
                    mask_pad=self.mask_pad,
                    snr_limit=self.snr_limit,
                    rms_limit=self.rms_limit,
                    fit_info=self.save_fit_info,
                )
                if self.save_fit_info:
                    result, fit_info = output
                    if fit_info is None:
                        log.info("No fit info to save; not creating output table.")
                    else:
                        table_filepath = self.make_output_path(ext=".ecsv", suffix="bfe_fit")
                        fit_info.write(table_filepath, format="ascii.ecsv", overwrite=True)
                        log.info(f"Saved fit info table in {table_filepath}")
                else:
                    result = output

        result.meta.cal_step.charge_migration = "COMPLETE"

        return result
