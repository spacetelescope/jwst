import logging

from stdatamodels.jwst import datamodels

from jwst.charge_migration import charge_migration
from jwst.stpipe import Step

log = logging.getLogger(__name__)

__all__ = ["ChargeMigrationStep"]


class ChargeMigrationStep(Step):
    """Set the CHARGELOSS flag for groups exhibiting significant charge migration."""

    class_alias = "charge_migration"

    spec = """
        signal_threshold = float(default=25000)
        skip = boolean(default=True)
    """  # noqa: E501

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

        # Retrieve the parameter value(s)
        signal_threshold = self.signal_threshold

        result = charge_migration.charge_migration(result, signal_threshold)
        result.meta.cal_step.charge_migration = "COMPLETE"

        return result
