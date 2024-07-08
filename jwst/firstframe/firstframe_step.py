import gc
from ..stpipe import Step
from . import firstframe_sub
from jwst.lib.basic_utils import use_datamodel, copy_datamodel


__all__ = ["FirstFrameStep"]


class FirstFrameStep(Step):
    """
    FirstFrameStep: This is a MIRI specific task.  If the number of groups
    is greater than 3, the DO_NOT_USE group data quality flag is added to
    first group.
    """

    class_alias = "firstframe"

    spec = """
    """

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model)

        result, input_model = copy_datamodel(input_model, self.parent)

        # check the data is MIRI data
        detector = result.meta.instrument.detector.upper()
        if detector[:3] == 'MIR':
            # Do the firstframe correction subtraction
            result = firstframe_sub.do_correction(result)
        else:
            self.log.warning('First Frame Correction is only for MIRI data')
            self.log.warning('First frame step will be skipped')
            result.meta.cal_step.firstframe = 'SKIPPED'

        gc.collect()
        return result
