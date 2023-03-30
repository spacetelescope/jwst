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
from ..rscd import rscd_step
from ..firstframe import firstframe_step
from ..lastframe import lastframe_step
from ..linearity import linearity_step
from ..dark_current import dark_current_step
from ..reset import reset_step
from ..persistence import persistence_step
from ..jump import jump_step
from ..undersampling_correction import  undersampling_correction_step
from ..ramp_fitting import ramp_fit_step
from ..gain_scale import gain_scale_step

__all__ = ['Detector1Pipeline']

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Detector1Pipeline(Pipeline):
    """
    Detector1Pipeline: Apply all calibration steps to raw JWST
    ramps to produce a 2-D slope product. Included steps are:
    group_scale, dq_init, saturation, ipc, superbias, refpix, rscd,
    lastframe, linearity, dark_current, persistence, jump detection,
    ramp_fit, and gain_scale.
    """

    class_alias = "calwebb_detector1"

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
                 'undersampling_correction': undersampling_correction_step.UndersamplingCorrectionStep,
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

        instrument = input.meta.instrument.name
        if instrument == 'MIRI':

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
            input = self.dark_current(input)
            input = self.refpix(input)

            # skip until MIRI team has figured out an algorithm
            # input = self.persistence(input)

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

            # skip persistence for NIRSpec
            if instrument != 'NIRSPEC':
                input = self.persistence(input)

            input = self.dark_current(input)

        # apply the jump step
        input = self.jump(input)

        # apply the undersampling_correction step
        input = self.undersampling_correction(input)

        # save the corrected ramp data, if requested
        if self.save_calibrated_ramp:
            self.save_model(input, 'ramp')

        # apply the ramp_fit step
        # This explicit test on self.ramp_fit.skip is a temporary workaround
        # to fix the problem that the ramp_fit step ordinarily returns two
        # objects, but when the step is skipped due to `skip = True`,
        # only the input is returned when the step is invoked.
        if self.ramp_fit.skip:
            input = self.ramp_fit(input)
            ints_model = None
        else:
            input, ints_model = self.ramp_fit(input)

        # apply the gain_scale step to the exposure-level product
        if input is not None:
            self.gain_scale.suffix = 'gain_scale'
            input = self.gain_scale(input)
        else:
            log.info("NoneType returned from ramp_fit.  Gain Scale step skipped.")

        # apply the gain scale step to the multi-integration product,
        # if it exists, and then save it
        if ints_model is not None:
            self.gain_scale.suffix = 'gain_scaleints'
            ints_model = self.gain_scale(ints_model)
            self.save_model(ints_model, 'rateints')

        # setup output_file for saving
        self.setup_output(input)

        log.info('... ending calwebb_detector1')

        return input

    def setup_output(self, input):
        if input is None:
            return None
        # Determine the proper file name suffix to use later
        if input.meta.cal_step.ramp_fit == 'COMPLETE':
            self.suffix = 'rate'
        else:
            self.suffix = 'ramp'
