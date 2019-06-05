#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import straylight

__all__ = ["StraylightStep"]


class StraylightStep (Step):
    """
    StraylightStep: Performs straylight correction image using a Mask file.
    """

    spec = """
         method = option('Nearest','ModShepard',default='ModShepard') #Algorithm method
         roi = float(default = 50.0) # Region of interest
         power = float(default = 1.0) # Power of weighting function

    """
    reference_file_types = ['regions']
    def process(self, input):


        # Open the input data model
        with datamodels.IFUImageModel(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if detector == 'MIRIFUSHORT':

                if self.method == 'Nearest':
                    # Use the Regions reference file set to 20% throughput threshhold
                    self.straylight_name = self.get_reference_file(input_model,
                                                                   'regions')
                    self.log.info('Using regions reference file %s',
                                  self.straylight_name)
                    # Check for a valid reference file
                    if self.straylight_name == 'N/A':
                        self.log.warning('No REGIONS reference file found')
                        self.log.warning('Straylight step will be skipped')
                        result = input_model.copy()
                        result.meta.cal_step.straylight = 'SKIPPED'
                        return result

                    allregions = datamodels.RegionsModel(self.straylight_name)
                    # Use 20% throughput array
                    regions=(allregions.regions)[2,:,:].copy()
                    self.log.info(' Using 20% throughput threshhold.')
                    self.log.info(' Using row-by-row approach.')
                    # Do the correction
                    result = straylight.correct_mrs(input_model, regions)
                    # Close the reference file and update the step status
                    allregions.close()
# ________________________________________________________________________________
                if self.method == 'ModShepard':
                    # Use the Regions reference file set to 20% throughput threshhold
                    self.straylight_name = self.get_reference_file(input_model,
                                                                   'regions')
                    self.log.info('Using regions reference file %s',
                                  self.straylight_name)
                    # Check for a valid reference file
                    if self.straylight_name == 'N/A':
                        self.log.warning('No REGIONS reference file found')
                        self.log.warning('Straylight step will be skipped')
                        result = input_model.copy()
                        result.meta.cal_step.straylight = 'SKIPPED'
                        return result

                    allregions = datamodels.RegionsModel(self.straylight_name)
                    # Use 20% throughput array
                    regions=(allregions.regions)[2,:,:].copy()
                    self.log.info(' Using 20% throughput threshhold.')
                    self.log.info(' Region of influence radius (pixels) %6.2f', self.roi)
                    self.log.info(' Modified Shepard weighting power %5.2f', self.power)
                    # Do the correction
                    result = straylight.correct_mrs_modshepard(input_model,
                                                               regions,
                                                               self.roi,
                                                               self.power)

                    # Close the reference file and update the step status
                    allregions.close()

# ________________________________________________________________________________
                result.meta.cal_step.straylight = 'COMPLETE'

            else:
                self.log.warning('Straylight correction not defined for detector %s',
                                  detector)
                self.log.warning('Straylight step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.straylight = 'SKIPPED'

        return result

class ErrorNoAssignWCS(Exception):
    pass
