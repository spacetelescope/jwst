#!/usr/bin/env python
import logging

from stdatamodels.jwst import datamodels

from jwst.stpipe import Pipeline

# step imports
from jwst.group_scale import group_scale_step
from jwst.dq_init import dq_init_step
from jwst.emicorr import emicorr_step
from jwst.saturation import saturation_step
from jwst.ipc import ipc_step
from jwst.superbias import superbias_step
from jwst.refpix import refpix_step
from jwst.rscd import rscd_step
from jwst.firstframe import firstframe_step
from jwst.lastframe import lastframe_step
from jwst.linearity import linearity_step
from jwst.dark_current import dark_current_step
from jwst.reset import reset_step
from jwst.persistence import persistence_step
from jwst.charge_migration import charge_migration_step
from jwst.jump import jump_step
from jwst.clean_flicker_noise import clean_flicker_noise_step
from jwst.ramp_fitting import ramp_fit_step
from jwst.gain_scale import gain_scale_step

__all__ = ["Detector1Pipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Detector1Pipeline(Pipeline):
    """
    Apply all calibration steps to raw JWST ramps to produce a 2-D slope product.

    Included steps are:
    group_scale, dq_init, saturation, ipc, superbias, refpix, rscd,
    lastframe, linearity, dark_current, persistence, jump detection,
    ramp_fit, and gain_scale.
    """

    class_alias = "calwebb_detector1"

    spec = """
        save_calibrated_ramp = boolean(default=False)
    """  # noqa: E501

    # Define aliases to steps
    step_defs = {
        "group_scale": group_scale_step.GroupScaleStep,
        "dq_init": dq_init_step.DQInitStep,
        "emicorr": emicorr_step.EmiCorrStep,
        "saturation": saturation_step.SaturationStep,
        "ipc": ipc_step.IPCStep,
        "superbias": superbias_step.SuperBiasStep,
        "refpix": refpix_step.RefPixStep,
        "rscd": rscd_step.RscdStep,
        "firstframe": firstframe_step.FirstFrameStep,
        "lastframe": lastframe_step.LastFrameStep,
        "linearity": linearity_step.LinearityStep,
        "dark_current": dark_current_step.DarkCurrentStep,
        "reset": reset_step.ResetStep,
        "persistence": persistence_step.PersistenceStep,
        "charge_migration": charge_migration_step.ChargeMigrationStep,
        "jump": jump_step.JumpStep,
        "clean_flicker_noise": clean_flicker_noise_step.CleanFlickerNoiseStep,
        "ramp_fit": ramp_fit_step.RampFitStep,
        "gain_scale": gain_scale_step.GainScaleStep,
    }

    # start the actual processing
    def process(self, input_data):
        """
        Run the Detector1 pipeline on the input data.

        Parameters
        ----------
        input_data : str or `~jwst.datamodels.RampModel`
            The input data to process.

        Returns
        -------
        `~jwst.datamodels.JwstDataModel`
            The calibrated data model.
        """
        log.info("Starting calwebb_detector1 ...")

        # open the input data as a RampModel
        input_data = datamodels.RampModel(input_data)

        # propagate output_dir to steps that might need it
        self.dark_current.output_dir = self.output_dir
        self.ramp_fit.output_dir = self.output_dir

        instrument = input_data.meta.instrument.name
        if instrument == "MIRI":
            # process MIRI exposures;
            # the steps are in a different order than NIR
            log.debug("Processing a MIRI exposure")

            input_data = self.group_scale.run(input_data)
            input_data = self.dq_init.run(input_data)
            input_data = self.emicorr.run(input_data)
            input_data = self.saturation.run(input_data)
            input_data = self.ipc.run(input_data)
            input_data = self.firstframe.run(input_data)
            input_data = self.lastframe.run(input_data)
            input_data = self.reset.run(input_data)
            input_data = self.linearity.run(input_data)
            input_data = self.rscd.run(input_data)
            input_data = self.dark_current.run(input_data)
            input_data = self.refpix.run(input_data)

            # skip until MIRI team has figured out an algorithm
            # input_data = self.persistence(input_data)

        else:
            # process Near-IR exposures
            log.debug("Processing a Near-IR exposure")

            input_data = self.group_scale.run(input_data)
            input_data = self.dq_init.run(input_data)
            input_data = self.saturation.run(input_data)
            input_data = self.ipc.run(input_data)
            input_data = self.superbias.run(input_data)
            input_data = self.refpix.run(input_data)
            input_data = self.linearity.run(input_data)

            # skip persistence for NIRSpec
            if instrument != "NIRSPEC":
                input_data = self.persistence.run(input_data)

            input_data = self.dark_current.run(input_data)

        # apply the charge_migration step
        input_data = self.charge_migration.run(input_data)

        # apply the jump step
        input_data = self.jump.run(input_data)

        # apply the clean_flicker_noise step
        input_data = self.clean_flicker_noise.run(input_data)

        # save the corrected ramp data, if requested
        if self.save_calibrated_ramp:
            self.save_model(input_data, "ramp")

        # apply the ramp_fit step
        # This explicit test on self.ramp_fit.skip is a temporary workaround
        # to fix the problem that the ramp_fit step ordinarily returns two
        # objects, but when the step is skipped due to `skip = True`,
        # only the input is returned when the step is invoked.
        if self.ramp_fit.skip:
            input_data = self.ramp_fit.run(input_data)
            ints_model = None
        else:
            input_data, ints_model = self.ramp_fit.run(input_data)

        # apply the gain_scale step to the exposure-level product
        if input_data is not None:
            self.gain_scale.suffix = "gain_scale"
            input_data = self.gain_scale.run(input_data)
        else:
            log.info("NoneType returned from ramp_fit.  Gain Scale step skipped.")

        # apply the gain scale step to the multi-integration product,
        # if it exists, and then save it
        if ints_model is not None:
            self.gain_scale.suffix = "gain_scaleints"
            ints_model = self.gain_scale.run(ints_model)
            self.save_model(ints_model, "rateints")

        # setup output_file for saving
        self.setup_output(input_data)

        log.info("... ending calwebb_detector1")

        return input_data

    def setup_output(self, input_data):
        """
        Set up the output file suffix based on which steps were run successfully.

        Parameters
        ----------
        input_data : `~jwst.datamodels.JwstDataModel`
            The output data product from the Detector1 pipeline
        """
        if input_data is None:
            return
        # Determine the proper file name suffix to use later
        if input_data.meta.cal_step.ramp_fit == "COMPLETE":
            self.suffix = "rate"
        else:
            self.suffix = "ramp"
