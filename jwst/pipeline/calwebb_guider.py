#!/usr/bin/env python
from ..stpipe import Pipeline
from .. import datamodels
import os

# step imports
from ..dq_init import dq_init_step
from ..flatfield import flat_field_step
from ..guider_cds import guider_cds_step

__version__ = "7.0.0"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class GuiderPipeline(Pipeline):
    """

    GuiderPipeline: For FGS observations, apply all calibration 
    steps to raw JWST ramps to produce a 2-D slope product. 
    Included steps are: dq_init, and flat_field.
    """

    spec = """
        save_calibrated_ramp = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {'dq_init': dq_init_step.DQInitStep,
                 'flat_field': flat_field_step.FlatFieldStep,
                 'guider_cds': guider_cds_step.GuiderCdsStep,
                 }


    # start the actual processing
    def process(self, input):

        log.info('Starting calwebb_guider ...')

        # open the input
        input = datamodels.open(input)

        # propagate output_dir to steps that might need it
        ## self.dark_current.output_dir = self.output_dir # dg - I think I don't want this
        ## self.ramp_fit.output_dir = self.output_dir # dg - I think I don't want this

        input = self.dq_init(input)  
        input = self.flat_field(input)
        input = self.guider_cds(input)

        # save the corrected ramp data, if requested
        ## if self.save_calibrated_ramp:
        ##    self.save_model(input, 'ramp')

        # setup output_file for saving
        self.setup_output(input)

        log.info('... ending calwebb_guider')

        return input


    def setup_output(self, input):

        # This routine doesn't actually save the final result to a file,
        # but just sets up the value of self.output_file appropriately.
        # The final data model is passed back up to the caller, which can be
        # either an interactive session or a command-line instance of stpipe.
        # If it's an interactive session, the data model is simply returned to
        # the user without saving to a file. If it's a command-line instance
        # of stpipe, stpipe will save the data model to a file using the name
        # given in self.output_file.

        # first determine the proper file name suffix to use later
        if input.meta.cal_step.ramp_fit == 'COMPLETE':
            suffix = 'rate'
        else:
            suffix = 'ramp'

        # Has an output file name already been set?
        if self.output_file is not None:

            # Check to see if the output_file name is the default set by
            # stpipe for command-line processing
            root, ext = os.path.splitext(self.output_file)
            if root[root.rfind('_') + 1:] == 'GuiderPipeline':

                # Remove the step name that stpipe appended to the file name,
                # as well as the original suffix on the input file name,
                # and create a new name with the appropriate output suffix
                root = root[:root.rfind('_')]
                self.output_file = root[:root.rfind('_') + 1] + suffix + ext

        # If no output name was set, take no action
