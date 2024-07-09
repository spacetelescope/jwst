#! /usr/bin/env python
import gc
import logging

from stdatamodels.jwst import datamodels
from ..stpipe import Step

from . import charge_migration
from jwst.lib.basic_utils import use_datamodel, copy_datamodel

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

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model, model_class=datamodels.RampModel)

        result, input_model = copy_datamodel(input_model, self.parent)

        if (result.data.shape[1] < 3):  # skip step if only 1 or 2 groups/integration
            log.info('Too few groups per integration; skipping charge_migration')

            result.meta.cal_step.charge_migration = 'SKIPPED'
            return result

        # Retrieve the parameter value(s)
        signal_threshold = self.signal_threshold

        result = charge_migration.charge_migration(result, signal_threshold)
        result.meta.cal_step.charge_migration = 'COMPLETE'

        gc.collect()
        return result
