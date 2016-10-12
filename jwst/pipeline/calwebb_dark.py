#!/usr/bin/env python
from ..stpipe import Pipeline
from .. import datamodels
import os

# step imports
from ..dq_init import dq_init_step
from ..saturation import saturation_step
from ..ipc import ipc_step
from ..superbias import superbias_step
from ..refpix import refpix_step
from ..rscd import rscd_step
from ..lastframe import lastframe_step
from ..linearity import linearity_step


__version__ = "0.7"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class DarkPipeline(Pipeline):
    """

    DarkPipeline: Apply detector-level calibration steps to raw JWST
    dark ramp to produce a corrected 4-D ramp product.
    Included steps are: dq_init, saturation, ipc, superbias, refpix,
    rscd, lastframe, and linearity.

    """

    # Define aliases to steps
    step_defs = {'dq_init': dq_init_step.DQInitStep,
                 'saturation': saturation_step.SaturationStep,
                 'ipc': ipc_step.IPCStep,
                 'superbias': superbias_step.SuperBiasStep,
                 'refpix': refpix_step.RefPixStep,
                 'rscd': rscd_step.RSCD_Step,
                 'lastframe': lastframe_step.LastFrameStep,
                 'linearity': linearity_step.LinearityStep,
                 }


    # start the actual processing
    def process(self, input):

        log.info('Starting calwebb_dark ...')

        # open the input
        input = datamodels.open(input)

        if input.meta.instrument.name == 'MIRI':

            # process MIRI exposures;
            # the steps are in a different order than NIR
            log.debug('Processing a MIRI exposure')

            input = self.dq_init(input)
            input = self.saturation(input)
            input = self.ipc(input)
            input = self.linearity(input)
            input = self.rscd(input)
            input = self.lastframe(input)

        else:

            # process Near-IR exposures
            log.debug('Processing a Near-IR exposure')

            input = self.dq_init(input)
            input = self.saturation(input)
            input = self.ipc(input)
            input = self.superbias(input)
            input = self.refpix(input)
            input = self.linearity(input)

        # setup output_file for saving
        self.setup_output(input)

        log.info('... ending calwebb_dark')

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

        # The desired suffix for corrected dark ramps is 'dark'
        suffix = 'dark'

        # Has an output file name already been set, either by the user
        # or by stpipe when the pipeline was instantiated?
        if self.output_file is not None:

            # Check to see if the output_file name is the default set by
            # stpipe for command-line processing
            root, ext = os.path.splitext(self.output_file)
            if root[root.rfind('_') + 1:] == 'DarkPipeline':

                # Remove the class name that stpipe appended to the file name,
                # as well as the original suffix on the input file name,
                # and create a new name with the appropriate output suffix
                root = root[:root.rfind('_')]
                self.output_file = root[:root.rfind('_') + 1] + suffix + ext

        # If no output name was set, take no action
