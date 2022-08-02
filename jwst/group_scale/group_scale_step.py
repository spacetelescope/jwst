from ..stpipe import Step
from .. import datamodels
from . import group_scale


__all__ = ["GroupScaleStep"]


class GroupScaleStep(Step):
    """
    GroupScaleStep: Rescales group data to account for on-board
    frame averaging that did not use FRMDIVSR = NFRAMES.
    All groups in the exposure are rescaled by FRMDIVSR/NFRAMES.
    """

    class_alias = "group_scale"

    def process(self, input):

        # Open the input data model
        with datamodels.RampModel(input) as input_model:

            # Always work on a copy
            result = input_model.copy()

            # Try to get values of NFRAMES and FRMDIVSR to see
            # if we need to do any rescaling
            nframes = result.meta.exposure.nframes
            frame_divisor = result.meta.exposure.frame_divisor

            # If we didn't find NFRAMES, we don't have enough info
            # to continue. Skip the step.
            if nframes is None:
                self.log.warning('NFRAMES value not found')
                self.log.warning('Step will be skipped')
                result.meta.cal_step.group_scale = 'SKIPPED'
                return result

            # If we didn't find FRMDIVSR, then check to see if NFRAMES
            # is a power of 2. If it is, rescaling isn't needed.
            if frame_divisor is None:
                if (nframes & (nframes - 1) == 0):
                    self.log.info('NFRAMES={} is a power of 2; correction not needed'.format(nframes))
                    self.log.info('Step will be skipped')
                    result.meta.cal_step.group_scale = 'SKIPPED'
                    return result

            # Compare NFRAMES and FRMDIVSR. If they're equal,
            # rescaling isn't needed.
            elif nframes == frame_divisor:
                self.log.info('NFRAMES and FRMDIVSR are equal; correction not needed')
                self.log.info('Step will be skipped')
                result.meta.cal_step.group_scale = 'SKIPPED'
                return result

            # Do the scaling
            group_scale.do_correction(result)

        return result
