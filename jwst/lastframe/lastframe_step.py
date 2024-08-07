from stdatamodels.jwst import datamodels
from ..stpipe import Step
from . import lastframe_sub

__all__ = ["LastFrameStep"]


class LastFrameStep(Step):
    """
    LastFrameStep: This is a MIRI specific task.  If the number of groups
    is greater than 2, the GROUP data quality flags for the final group will
    be set to DO_NOT_USE.
    """

    class_alias = "lastframe"

    spec = """
    """

    def process(self, step_input):

        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector

            if detector[:3] != 'MIR':
                self.log.warning('Last Frame Correction is only for MIRI data')
                self.log.warning('Last frame step will be skipped')
                input_model.meta.cal_step.lastframe = 'SKIPPED'
                return input_model

            # Work on a copy
            result = input_model.copy()

            # Do the lastframe correction subtraction
            result = lastframe_sub.do_correction(result)

        return result
