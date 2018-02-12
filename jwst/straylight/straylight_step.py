#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import straylight
from ..datamodels import RegionsModel

class StraylightStep (Step):
    """
    StraylightStep: Performs straylight correction image using a Mask file.
    """

    spec = """
         method = option('Nearest','ModShepard',default='ModShepard')
         roi = float(default = 50.0)
         power = float(default = 1.0)

    """
    reference_file_types = ['regions','straymask']

    def process(self, input):


        # Open the input data model
        with datamodels.IFUImageModel(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if detector == 'MIRIFUSHORT':

                if self.method == 'Nearest': 
                # Get the name of the straylight reference file
                    self.straylight_name = self.get_reference_file(input_model,
                                                                   'straymask')
                    self.log.info('Using straylight reference file %s',
                                  self.straylight_name)
                # Check for a valid reference file
                    if self.straylight_name == 'N/A':
                        self.log.warning('No STRAYLIGHT reference file found')
                        self.log.warning('Straylight step will be skipped')
                        result = input_model.copy()
                        result.meta.cal_step.straylight = 'SKIPPED'
                        return result

                # Open the straylight mask ref file data model
                    straylight_model = datamodels.StrayLightModel(self.straylight_name)
                    result = straylight.correct_MRS(input_model, straylight_model)
                # Close the reference file and update the step status
                    straylight_model.close()
# ________________________________________________________________________________
                if self.method == 'ModShepard': 
                    self.regions_name = self.get_reference_file(input_model,
                                                                'regions')
                    self.log.info('Using Regions reference file %s',
                                  self.regions_name)
                # Check for a valid reference file
                    if self.regions_name == 'N/A':
                        self.log.warning('No STRAYLIGHT/REGIONS reference file found')
                        self.log.warning('Straylight step will be skipped')
                        result = input_model.copy()
                        result.meta.cal_step.straylight = 'SKIPPED'
                        return result

                # Open the straylight mask ref file data model
                    region_model = datamodels.RegionsModel(self.regions_name)

                # Do the correction
                    result = straylight.correct_MRS_ModShepard(input_model, 
                                                               region_model,
                                                               self.roi,
                                                               self.power)
# ________________________________________________________________________________
                result.meta.cal_step.straylight = 'COMPLETE'

            else:
                self.log.warning('Straylight correction not defined for detector %s',
                                  detector)
                self.log.warning('Straylight step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.straylight = 'SKIPPED'

        return result

