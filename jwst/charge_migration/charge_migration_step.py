#! /usr/bin/env python
import logging

from stdatamodels.jwst import datamodels
from ..stpipe import Step

from . import charge_migration

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ChargeMigrationStep"]


class ChargeMigrationStep(Step):
    """
    This Step sets the CHARGELOSS flag for groups exhibiting significant
    charge migration.
    """
    class_alias = "charge_migration"

    spec = """
        signal_threshold = float(default=25000)
        skip = boolean(default=True)
    """

    def process(self, step_input):

        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:

            if (input_model.data.shape[1] < 3):  # skip step if only 1 or 2 groups/integration
                log.info('Too few groups per integration; skipping charge_migration')

                input_model.meta.cal_step.charge_migration = 'SKIPPED'
                return input_model

            # Work on a copy
            result = input_model.copy()

            # Retrieve the parameter value(s)
            signal_threshold = self.signal_threshold

            result = charge_migration.charge_migration(result, signal_threshold)
            result.meta.cal_step.charge_migration = 'COMPLETE'

        return result
