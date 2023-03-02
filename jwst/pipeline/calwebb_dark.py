#!/usr/bin/env python
import logging

from stdatamodels.jwst import datamodels

from ..stpipe import Pipeline

# step imports
from ..group_scale import group_scale_step
from ..dq_init import dq_init_step
from ..saturation import saturation_step
from ..ipc import ipc_step
from ..superbias import superbias_step
from ..refpix import refpix_step
from ..reset import reset_step
from ..rscd import rscd_step
from ..firstframe import firstframe_step
from ..lastframe import lastframe_step
from ..linearity import linearity_step

__all__ = ['DarkPipeline']

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class DarkPipeline(Pipeline):
    """

    DarkPipeline: Apply detector-level calibration steps to raw JWST
    dark ramp to produce a corrected 4-D ramp product.
    Included steps are: group_scale, dq_init, saturation, ipc,
    superbias, refpix, rscd, lastframe, and linearity.

    """

    class_alias = "calwebb_dark"

    # Define aliases to steps
    step_defs = {'group_scale': group_scale_step.GroupScaleStep,
                 'dq_init': dq_init_step.DQInitStep,
                 'saturation': saturation_step.SaturationStep,
                 'ipc': ipc_step.IPCStep,
                 'superbias': superbias_step.SuperBiasStep,
                 'refpix': refpix_step.RefPixStep,
                 'reset': reset_step.ResetStep,
                 'rscd': rscd_step.RscdStep,
                 'firstframe': firstframe_step.FirstFrameStep,
                 'lastframe': lastframe_step.LastFrameStep,
                 'linearity': linearity_step.LinearityStep,
                 }

    # start the actual processing
    def process(self, input):

        log.info('Starting calwebb_dark ...')

        # open the input
        input = datamodels.RampModel(input)

        if input.meta.instrument.name == 'MIRI':

            # process MIRI exposures;
            # the steps are in a different order than NIR
            log.debug('Processing a MIRI exposure')

            input = self.group_scale(input)
            input = self.dq_init(input)
            input = self.saturation(input)
            input = self.ipc(input)
            input = self.firstframe(input)
            input = self.lastframe(input)
            input = self.reset(input)
            input = self.linearity(input)
            input = self.rscd(input)

        else:

            # process Near-IR exposures
            log.debug('Processing a Near-IR exposure')

            input = self.group_scale(input)
            input = self.dq_init(input)
            input = self.saturation(input)
            input = self.ipc(input)
            input = self.superbias(input)
            input = self.refpix(input)
            input = self.linearity(input)

        log.info('... ending calwebb_dark')

        return input
