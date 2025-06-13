#! /usr/bin/env python


from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import guider_cds

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["GuiderCdsStep"]


class GuiderCdsStep(Step):
    """Calculate the count rate for each pixel for FGS modes."""

    class_alias = "guider_cds"

    def process(self, input_data):
        """
        Execute the step.

        Parameters
        ----------
        input_data : datamodel, str
            The input GuiderRawModel or filename containing
            a GuiderRawModel.

        Returns
        -------
        JWST GuiderCalModel or GuiderRawModel
            This will be a GuiderRawModel if the step was skipped; otherwise,
            it will be a GuiderCalModel containing calibrated count rates.
        """
        with datamodels.GuiderRawModel(input_data) as input_model:
            # Get the gain reference file
            gain_filename = self.get_reference_file(input_model, "gain")
            if gain_filename == "N/A":
                self.log.warning("No GAIN reference file found!")
                self.log.warning("guider_cds step will be skipped.")
                input_model.meta.cal_step.guider_cds = "SKIPPED"
                return input_model

            self.log.info("Using GAIN reference file: %s", gain_filename)
            gain_model = datamodels.GainModel(gain_filename)

            # Get the readnoise reference file
            readnoise_filename = self.get_reference_file(input_model, "readnoise")
            if readnoise_filename == "N/A":
                self.log.warning("No READNOISE reference file found!")
                self.log.warning("guider_cds step will be skipped.")
                input_model.meta.cal_step.guider_cds = "SKIPPED"
                return input_model

            self.log.info("Using READNOISE reference file: %s", readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)

            out_model = guider_cds.guider_cds(input_model, gain_model, readnoise_model)

        out_model.meta.cal_step.guider_cds = "COMPLETE"

        return out_model
