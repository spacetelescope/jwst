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
from jwst.reset import reset_step
from jwst.rscd import rscd_step
from jwst.firstframe import firstframe_step
from jwst.lastframe import lastframe_step
from jwst.linearity import linearity_step

__all__ = ["DarkPipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class DarkPipeline(Pipeline):
    """
    Produce a 4-D corrected ramp product from a raw JWST dark ramp.

    Included steps are: group_scale, dq_init, saturation, ipc,
    superbias, refpix, rscd, lastframe, and linearity.
    """

    class_alias = "calwebb_dark"

    # Define aliases to steps
    step_defs = {
        "group_scale": group_scale_step.GroupScaleStep,
        "dq_init": dq_init_step.DQInitStep,
        "emicorr": emicorr_step.EmiCorrStep,
        "saturation": saturation_step.SaturationStep,
        "ipc": ipc_step.IPCStep,
        "superbias": superbias_step.SuperBiasStep,
        "refpix": refpix_step.RefPixStep,
        "reset": reset_step.ResetStep,
        "rscd": rscd_step.RscdStep,
        "firstframe": firstframe_step.FirstFrameStep,
        "lastframe": lastframe_step.LastFrameStep,
        "linearity": linearity_step.LinearityStep,
    }

    # start the actual processing
    def process(self, input_data):
        """
        Run the calwebb_dark pipeline on the input data.

        Parameters
        ----------
        input_data : str or `~jwst.datamodels.RampModel`
            The input data to process. If a string, it is assumed
            to be a filename pointing to a RampModel.

        Returns
        -------
        RampModel
            The 4-D corrected ramp product.
        """
        log.info("Starting calwebb_dark ...")

        # open the input
        input_data = datamodels.RampModel(input_data)

        if input_data.meta.instrument.name == "MIRI":
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

        log.info("... ending calwebb_dark")

        return input_data
