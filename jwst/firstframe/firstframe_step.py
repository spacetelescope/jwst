import logging

from stdatamodels.jwst import datamodels

from jwst.firstframe import firstframe_sub
from jwst.stpipe import Step

__all__ = ["FirstFrameStep"]

log = logging.getLogger(__name__)


class FirstFrameStep(Step):
    """
    Set data quality flags for the first group in MIRI ramps.

    A MIRI specific task.  If the number of groups is > than 3,
    the DO_NOT_USE group data quality flag is added to first group.
    """

    class_alias = "firstframe"

    spec = """
        bright_use_group1 = boolean(default=False) # do not flag group1 if group3 is saturated
    """  # noqa: E501

    def process(self, step_input):
        """
        For MIRI data with more than 3 groups, set first group dq to DO_NOT_USE.

        Parameters
        ----------
        step_input : DataModel
            Input datamodel to be corrected

        Returns
        -------
        output_model : DataModel
            Firstframe corrected datamodel
        """
        # Open the input data model
        with datamodels.open(step_input) as input_model:
            # check the data is MIRI data
            detector = input_model.meta.instrument.detector.upper()
            if detector[:3] != "MIR":
                log.warning("First Frame Correction is only for MIRI data")
                log.warning("First frame step will be skipped")
                input_model.meta.cal_step.firstframe = "SKIPPED"
                return input_model

            # Cork on a copy
            result = input_model.copy()

            # Do the firstframe correction subtraction
            result = firstframe_sub.do_correction(result, bright_use_group1=self.bright_use_group1)

        return result
