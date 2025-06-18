#!/usr/bin/env python
import logging

from stdatamodels.jwst import datamodels

from jwst.stpipe import Pipeline

# step imports
from jwst.dq_init import dq_init_step
from jwst.flatfield import flat_field_step
from jwst.guider_cds import guider_cds_step

__all__ = ["GuiderPipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class GuiderPipeline(Pipeline):
    """
    Apply all calibration steps to raw JWST FGS observation ramps to produce a 3-D slope product.

    Included steps are: dq_init, guider_cds, and flat_field.
    """

    class_alias = "calwebb_guider"

    # Define aliases to steps
    step_defs = {
        "dq_init": dq_init_step.DQInitStep,
        "guider_cds": guider_cds_step.GuiderCdsStep,
        "flat_field": flat_field_step.FlatFieldStep,
    }

    # Start the processing
    def process(self, input_data):
        """
        Run the Guider1 pipeline on the input data.

        Parameters
        ----------
        input_data : str, GuiderRawModel, or GuiderCalModel
            The input data to process. If a string, it is assumed to be the
            name of a file to load. If a GuiderCalModel, the dq_init and guider_cds
            steps should be set to skip.

        Returns
        -------
        GuiderCalModel
            The calibrated data model.
        """
        # Set the output product type
        self.suffix = "cal"

        log.info("Starting calwebb_guider ...")

        # Open the input:
        # If the first two steps are set to be skipped, assume
        # they've been run before and open the input as a Cal
        # model, appropriate for input to flat_field
        if self.dq_init.skip and self.guider_cds.skip:
            log.info(
                "dq_init and guider_cds are set to skip; assume they"
                " were run before and load data as GuiderCalModel"
            )
            input_data = datamodels.GuiderCalModel(input_data)
        else:
            input_data = datamodels.GuiderRawModel(input_data)

        # Apply the steps
        input_data = self.dq_init.run(input_data)
        input_data = self.guider_cds.run(input_data)
        input_data = self.flat_field.run(input_data)

        log.info("... ending calwebb_guider")

        return input_data
