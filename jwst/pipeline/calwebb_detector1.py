#!/usr/bin/env python
import logging
from ..stpipe import Pipeline
from .. import datamodels

# step imports
from ..group_scale import group_scale_step
from ..dq_init import dq_init_step
from ..saturation import saturation_step
from ..ipc import ipc_step
from ..superbias import superbias_step
from ..refpix import refpix_step
from ..rscd import rscd_step
from ..firstframe import firstframe_step
from ..lastframe import lastframe_step
from ..linearity import linearity_step
from ..dark_current import dark_current_step
from ..reset import reset_step
from ..persistence import persistence_step
from ..jump import jump_step
from ..ramp_fitting import ramp_fit_step
from ..gain_scale import gain_scale_step

__all__ = ['Detector1Pipeline']

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class Detector1Pipeline(Pipeline):
    """
    Detector1Pipeline: Apply all calibration steps to raw JWST
    ramps to produce a 2-D slope product. Included steps are:
    group_scale, dq_init, saturation, ipc, superbias, refpix, rscd,
    lastframe, linearity, dark_current, persistence, jump detection,
    ramp_fit, and gain_scale.
    """

    spec = """
        save_calibrated_ramp = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {'group_scale': group_scale_step.GroupScaleStep,
                 'dq_init': dq_init_step.DQInitStep,
                 'saturation': saturation_step.SaturationStep,
                 'ipc': ipc_step.IPCStep,
                 'superbias': superbias_step.SuperBiasStep,
                 'refpix': refpix_step.RefPixStep,
                 'rscd': rscd_step.RscdStep,
                 'firstframe': firstframe_step.FirstFrameStep,
                 'lastframe': lastframe_step.LastFrameStep,
                 'linearity': linearity_step.LinearityStep,
                 'dark_current': dark_current_step.DarkCurrentStep,
                 'reset': reset_step.ResetStep,
                 'persistence': persistence_step.PersistenceStep,
                 'jump': jump_step.JumpStep,
                 'ramp_fit': ramp_fit_step.RampFitStep,
                 'gain_scale': gain_scale_step.GainScaleStep,
                 }

    # start the actual processing
    def process(self, input):

        log.info('Starting calwebb_detector1 ...')

        # open the input as a RampModel
        input = datamodels.RampModel(input)

        # propagate output_dir to steps that might need it
        self.dark_current.output_dir = self.output_dir
        self.ramp_fit.output_dir = self.output_dir

        if input.meta.instrument.name == 'MIRI':

            # process MIRI exposures;
            # the steps are in a different order than NIR
            log.debug('Processing a MIRI exposure')

            result = self.group_scale(input)
            result = self.dq_init(result)
            result = self.saturation(result)
            result = self.ipc(result)
            result = self.firstframe(result)
            result = self.lastframe(result)
            result = self.reset(result)
            result = self.linearity(result)
            result = self.rscd(result)
            result = self.dark_current(result)
            result = self.refpix(result)

            # skip until MIRI team has figured out an algorithm
            #result = self.persistence(result)

        else:

            # process Near-IR exposures
            log.debug('Processing a Near-IR exposure')

            result = self.group_scale(input)
            result = self.dq_init(result)
            result = self.saturation(result)
            result = self.ipc(result)
            result = self.superbias(result)
            result = self.refpix(result)
            result = self.linearity(result)

            # skip persistence for NIRSpec
            if result.meta.instrument.name != 'NIRSPEC':
                result = self.persistence(result)

            result = self.dark_current(result)

        # apply the jump step
        result = self.jump(result)

        # save the corrected ramp data, if requested
        if self.save_calibrated_ramp:
            result.meta.filetype = 'calibrated ramp'
            self.save_model(result, 'ramp')

        # apply the ramp_fit step
        # This explicit test on self.ramp_fit.skip is a temporary workaround
        # to fix the problem that the ramp_fit step ordinarily returns two
        # objects, but when the step is skipped due to `skip = True` in a
        # cfg file, only the input is returned when the step is invoked.
        if self.ramp_fit.skip:
            result = self.ramp_fit(result)
            ints_model = None
        else:
            result, ints_model = self.ramp_fit(result)

        # apply the gain_scale step to the exposure-level product
        self.gain_scale.suffix = 'gain_scale'
        result = self.gain_scale(result)

        # apply the gain scale step to the multi-integration product,
        # if it exists, and then save it
        if ints_model is not None:
            self.gain_scale.suffix = 'gain_scaleints'
            ints_model = self.gain_scale(ints_model)
            ints_model.meta.filetype = 'countrate'
            self.save_model(ints_model, 'rateints')

        # setup output_file for saving
        self.setup_output(result)
        result.meta.filetype = 'countrate'

        log.info('... ending calwebb_detector1')

        return result

    def setup_output(self, input):
        # Determine the proper file name suffix to use later
        if input.meta.cal_step.ramp_fit == 'COMPLETE':
            self.suffix = 'rate'
        else:
            self.suffix = 'ramp'
