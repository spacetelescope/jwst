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
        signal_threshold = float(default=10000)  # Threshold value for flagging charge loss.
        mask_pad = integer(default=2)  # Number of pixels to mask on either side of a pixel with charge loss.
        apply_correction = boolean(default=False)  # Fit and correct for charge migration.
        skip = boolean(default=True)  # By default, skip the step.
    """  # noqa: E501

    reference_file_types = ["readnoise"]

    def process(self, step_input):
        """
        Detect jumps based on signal threshold.

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

        if self.apply_correction:
            readnoise_filename = self.get_reference_file(result, "readnoise")
            log.info("Using READNOISE reference file: %s", readnoise_filename)

            with datamodels.ReadnoiseModel(readnoise_filename) as rnoise:
                result = correct_charge_migration(
                    result,
                    rnoise,
                    signal_threshold=self.signal_threshold,
                    mask_pad=self.mask_pad,
                )

        else:
            result = flag_charge_migration(result, self.signal_threshold)

        result.meta.cal_step.charge_migration = "COMPLETE"

        return result
