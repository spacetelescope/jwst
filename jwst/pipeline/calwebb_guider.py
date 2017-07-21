#!/usr/bin/env python
from ..stpipe import Pipeline
from .. import datamodels
import os

# step imports
from ..dq_init import dq_init_step
from ..flatfield import flat_field_step
from ..guider_cds import guider_cds_step

__version__ = "0.7.1"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class GuiderPipeline(Pipeline):
    """
    GuiderPipeline: For FGS observations, apply all calibration 
    steps to raw JWST ramps to produce a 3-D slope product. 
    Included steps are: dq_init, guider_cds, and flat_field.
    """

    # Define aliases to steps
    step_defs = {'dq_init': dq_init_step.DQInitStep,
                 'guider_cds': guider_cds_step.GuiderCdsStep,
                 'flat_field': flat_field_step.FlatFieldStep,
                 }

    # start the actual processing
    def process(self, input):

        log.info('Starting calwebb_guider ...')

        # open the input
        input = datamodels.GuiderRawModel(input) 
        input = self.dq_init(input)  
        input = self.guider_cds(input)
        input = self.flat_field(input)

        log.info('... ending calwebb_guider')

        return input
