#! /usr/bin/env python
import logging

from ..stpipe import Step

from stdatamodels.jwst import datamodels

from . import undersampling_correction

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["UndersamplingCorrectionStep"]


class UndersamplingCorrectionStep(Step):
    """
    This Step sets the UNDERSAMP flag for groups exhibiting significant
    charge migration.
    """
    class_alias = "undersampling_correction"

    spec = """
        signal_threshold = float(default=30000)
        skip = boolean(default=True)
    """

    def process(self, input):

        # Open the input data model
        with datamodels.RampModel(input) as input_model:
            if (input_model.data.shape[1] < 3):  # skip step if only 1 or 2 groups/integration
                log.info('Too few groups per integration; skipping undersampling_correction')
                
                result = input_model
                result.meta.cal_step.undersampling_correction = 'SKIPPED'

                return result

            # Retrieve the parameter value(s)
            signal_threshold = self.signal_threshold

            result = undersampling_correction.undersampling_correction(input_model, signal_threshold)
            result.meta.cal_step.undersampling_correction = 'COMPLETE'

        return result
