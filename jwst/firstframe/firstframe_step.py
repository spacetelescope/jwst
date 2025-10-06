import logging

from stdatamodels.jwst import datamodels

from jwst.firstframe import firstframe_sub
from jwst.stpipe import Step

__all__ = ["FirstFrameStep"]

log = logging.getLogger(__name__)


class FirstFrameStep(Step):
    """
    Set data quality flags for the first group in MIRI ramps.

    If the number of groups is > than 3, the DO_NOT_USE group data
    quality flag is added to first group.
    """

    class_alias = "firstframe"

    spec = """
        bright_use_group1 = boolean(default=True) # do not flag group1 if group3 is saturated
    """  # noqa: E501

    def process(self, step_input):
        """
        For MIRI data with more than 3 groups, set first group dq to DO_NOT_USE.

        Parameters
        ----------
        step_input : str or `~stdatamodels.jwst.datamodels.RampModel`
            Input filename or datamodel to be corrected.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.RampModel`
            First frame corrected datamodel.
        """
        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        # Check the data is MIRI data
        detector = result.meta.instrument.detector.upper()
        if detector[:3] != "MIR":
            log.warning("First Frame Correction is only for MIRI data")
            log.warning("First frame step will be skipped")
            result.meta.cal_step.firstframe = "SKIPPED"
            return result

        # Do the firstframe correction subtraction
        result = firstframe_sub.do_correction(result, bright_use_group1=self.bright_use_group1)

        return result
