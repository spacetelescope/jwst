#!/usr/bin/env python
from jwst.stpipe import Pipeline
from .. import datamodels
import os

# calwebb IMAGE2 step imports
from jwst.assign_wcs import assign_wcs_step
from jwst.flatfield import flat_field_step
from jwst.persistence import persistence_step
from jwst.emission import emission_step
from jwst.photom import photom_step


__version__ = "2.1"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class Image2Pipeline(Pipeline):
    """

    CalWebbImage2: Processes JWST imaging-mode slope images from 
                   Level-2a to Level-2b.

    Included steps are:
    assign_wcs, flat_field, persistence, emission, and photom.

    """

    # Define alias to steps
    step_defs = {'assign_wcs': assign_wcs_step.AssignWcsStep,
                 'flat_field': flat_field_step.FlatFieldStep,
                 'persistence': persistence_step.PersistenceStep,
                 'emission': emission_step.EmissionStep,
                 'photom': photom_step.PhotomStep,
                 }

    def process(self, input):

        log.info('Starting calwebb_image2 ...')

        # work on slope images
        input = self.assign_wcs(input)
        input = self.flat_field(input)
        input = self.persistence(input)
        input = self.emission(input)
        input = self.photom(input)

        # setup output_file for saving
        self.setup_output(input)

        log.info('... ending calwebb_image2')

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
        if isinstance(input, datamodels.CubeModel):
            suffix = 'calints'
        else:
            suffix = 'cal'

        # Has an output file name already been set?
        if self.output_file is not None:

            # Check to see if the output_file name is the default set by
            # stpipe for command-line processing
            root, ext = os.path.splitext(self.output_file)
            if root[root.rfind('_') + 1:] == 'Image2Pipeline':

                # Remove the step name that stpipe appended to the file name,
                # as well as the original suffix on the input file name,
                # and create a new name with the appropriate output suffix
                root = root[:root.rfind('_')]
                self.output_file = root[:root.rfind('_') + 1] + suffix + ext

        # If no output name was set, take no action
