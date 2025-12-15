import logging

from stdatamodels.jwst import datamodels

from jwst.firstframe import firstframe_sub
from jwst.stpipe import Step

__all__ = ["FirstFrameStep"]

log = logging.getLogger(__name__)


class FirstFrameStep(Step):
    """
    Set data quality flags for the first group in MIRI ramps.

    .. deprecated:: 1.21.0
        The `FirstFrameStep` has been deprecated and will be removed
        in a future release. Flagging the first grouops  has been added to the RSCD step.
    """

    class_alias = "firstframe"

    spec = """
    skip = booelan(default=True) # Do not run this step. 
    """  # noqa: E501

    def __init__(self, *args, **kwargs):
        deprecation_message = (
            "'FiristFrameStep' has been deprecated since 1.21.0 and "
            "will be removed in a future release. Flagging the first group  has been"
            " added to the RSCD step.  "
        )
        warnings.warn(deprecation_message, DeprecationWarning, stacklevel=2)
        log.warning(deprecation_message)
        super().__init__(*args, **kwargs)

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
