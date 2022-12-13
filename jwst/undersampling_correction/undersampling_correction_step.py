#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import undersampling_correction

__all__ = ["UndersamplingCorrectionStep"]


class UndersamplingCorrectionStep(Step):
    """
    This Step sets DO_NOT_USE flag for groups exhibiting significant
    charge migration.
    """
    class_alias = "undersampling_correction"

    spec = """
        signal_threshold = float(default=30000)
    """
    def process(self, input):

        # Open the input data model
        with datamodels.RampModel(input) as input_model:

            # Retrieve the parameter value(s)
            signal_threshold = self.signal_threshold

            result = undersampling_correction.undersampling_correction(input_model, signal_threshold)
            result.meta.cal_step.undersampling_correction = 'COMPLETE'

        return result
