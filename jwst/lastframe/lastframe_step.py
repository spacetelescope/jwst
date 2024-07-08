import gc
from stdatamodels.jwst import datamodels
from ..stpipe import Step
from . import lastframe_sub
from jwst.lib.basic_utils import use_datamodel, copy_datamodel

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

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model, model_class=datamodels.RampModel)

        result, input_model = copy_datamodel(input_model, self.parent)

        # check the data is MIRI data
        detector = result.meta.instrument.detector
        if detector[:3] == 'MIR':
            # Do the lastframe correction subtraction
            result = lastframe_sub.do_correction(result)
        else:
            self.log.warning('Last Frame Correction is only for MIRI data')
            self.log.warning('Last frame step will be skipped')
            result.meta.cal_step.lastframe = 'SKIPPED'

        gc.collect()
        return result
