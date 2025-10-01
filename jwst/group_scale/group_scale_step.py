import logging
import warnings

from stdatamodels.exceptions import ValidationWarning
from stdatamodels.jwst import datamodels

from jwst.group_scale import group_scale
from jwst.stpipe import Step

__all__ = ["GroupScaleStep"]

log = logging.getLogger(__name__)


class GroupScaleStep(Step):
    """Rescale group data to account for on-board frame averaging."""

    class_alias = "group_scale"

    spec = """
    """  # noqa: E501

    def process(self, step_input):
        """
        Perform group scale step.

        Rescales group data to account for on-board
        frame averaging that did not use FRMDIVSR = NFRAMES.
        All groups in the exposure are rescaled by FRMDIVSR/NFRAMES.

        Parameters
        ----------
        step_input : str or `~stdatamodels.jwst.datamodels.RampModel`
            Input file name or data model on which to perform group scale step.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            Output data model on which the group scale step has been performed.
        """
        # Open the input data model as a RampModel, catching a specific
        # expected warning for uncal files
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "error",
                    message=r"(?s:.*)Array datatype .* not compatible",
                    category=ValidationWarning,
                )
                result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)
        except ValidationWarning as err:
            log.error(err)

            # Inform the user and raise a clearer error message
            msg = (
                "Input data model does not have float-type data. "
                "The file should be opened as a RampModel before calling the step."
            )
            log.error(msg)
            raise TypeError(msg) from None

        # Try to get values of NFRAMES and FRMDIVSR to see
        # if we need to do any rescaling
        nframes = result.meta.exposure.nframes
        frame_divisor = result.meta.exposure.frame_divisor

        # If we didn't find NFRAMES, we don't have enough info
        # to continue. Skip the step.
        if nframes is None:
            log.warning("NFRAMES value not found")
            log.warning("Step will be skipped")
            result.meta.cal_step.group_scale = "SKIPPED"
            return result

        # If we didn't find FRMDIVSR, then check to see if NFRAMES
        # is a power of 2. If it is, rescaling isn't needed.
        if frame_divisor is None:
            if nframes & (nframes - 1) == 0:
                log.info(f"NFRAMES={nframes} is a power of 2; correction not needed")
                log.info("Step will be skipped")
                result.meta.cal_step.group_scale = "SKIPPED"
                return result

        # Compare NFRAMES and FRMDIVSR. If they're equal,
        # rescaling isn't needed.
        elif nframes == frame_divisor:
            log.info("NFRAMES and FRMDIVSR are equal; correction not needed")
            log.info("Step will be skipped")
            result.meta.cal_step.group_scale = "SKIPPED"
            return result

        # Do the scaling
        group_scale.do_correction(result)

        return result
