from stdatamodels.jwst import datamodels

from jwst.group_scale import group_scale
from jwst.stpipe import Step

__all__ = ["GroupScaleStep"]


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
        step_input : datamodel
            Input data model on which to perform group scale step.

        Returns
        -------
        result : datamodel
            Output data model on which the group scale step has been performed.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # Work on a copy
            result = input_model.copy()

            # Try to get values of NFRAMES and FRMDIVSR to see
            # if we need to do any rescaling
            nframes = result.meta.exposure.nframes
            frame_divisor = result.meta.exposure.frame_divisor

            # If we didn't find NFRAMES, we don't have enough info
            # to continue. Skip the step.
            if nframes is None:
                self.log.warning("NFRAMES value not found")
                self.log.warning("Step will be skipped")
                result.meta.cal_step.group_scale = "SKIPPED"
                return result

            # If we didn't find FRMDIVSR, then check to see if NFRAMES
            # is a power of 2. If it is, rescaling isn't needed.
            if frame_divisor is None:
                if nframes & (nframes - 1) == 0:
                    self.log.info(f"NFRAMES={nframes} is a power of 2; correction not needed")
                    self.log.info("Step will be skipped")
                    result.meta.cal_step.group_scale = "SKIPPED"
                    return result

            # Compare NFRAMES and FRMDIVSR. If they're equal,
            # rescaling isn't needed.
            elif nframes == frame_divisor:
                self.log.info("NFRAMES and FRMDIVSR are equal; correction not needed")
                self.log.info("Step will be skipped")
                result.meta.cal_step.group_scale = "SKIPPED"
                return result

            # Do the scaling
            group_scale.do_correction(result)

            # Cleanup
            del input_model

        return result
