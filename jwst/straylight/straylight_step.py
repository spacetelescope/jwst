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
         roi = integer(2, 1024, default = 50) # Region of interest given as even integer
         power = float(0.1, 5, default = 1.0) # Power of weighting function

    """

    reference_file_types = ['regions']

    def process(self, input):
        # Open the input data model
        with datamodels.IFUImageModel(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector

            if detector == 'MIRIFUSHORT':
                # If Modified Shepard  test  input parameters for weighting
                if self.method == 'ModShepard':
                    # reasonable power varible defined as: 0.1 < power < 5
                    if self.power < 0.1 or self.power > 5:
                        self.log.warning("The kernel power parameter is outside the reasonable range of"
                                  " 0.1 to 5. It is set to {}".format(self.power))
                        self.log.warning('Straylight step will be skipped')
                        result = input_model.copy()
                        result.meta.cal_step.straylight = 'SKIPPED'
                        return result

                    if self.roi < 2 or self.roi > 1024:
                        self.log.warning("The kernel roi parameter is outside the reasonable range of"
                                  " 2 to 1024. It is set to {} ".format(self.roi))
                        self.log.warning('Straylight step will be skipped')
                        result = input_model.copy()
                        result.meta.cal_step.straylight = 'SKIPPED'
                        return result

                    # test that ROI is an even number
                    # so that the kernel will be odd rows,columns in size
                    test = self.roi % 2
                    if test != 0:
                        self.log.info("The kernel roi parameter is odd value {}"
                                      "must be even. Adding 1 to size ".format(self.roi))
                        self.roi = self.roi + 1
                        if self.roi > 1024:
                            self.roi = self.roi - 2

                # Check for a valid reference file
                # Use the Regions reference file set to 20% throughput threshhold
                self.straylight_name = self.get_reference_file(input_model,
                                                               'regions')
                self.log.info('Using regions reference file %s',
                              self.straylight_name)

                if self.straylight_name == 'N/A':
                    self.log.warning('No REGIONS reference file found')
                    self.log.warning('Straylight step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.straylight = 'SKIPPED'
                    return result
                allregions = datamodels.RegionsModel(self.straylight_name)
                # Use 20% throughput array
                regions = (allregions.regions)[2, :, :].copy()
                self.log.info(' Using 20% throughput threshhold.')
                allregions.close()

                if self.method == 'Nearest':
                    # Do the correction
                    self.log.info(' Using row-by-row approach.')
                    result = straylight.correct_mrs(input_model, regions)
                elif self.method == 'ModShepard':
                    # Do the correction
                    self.log.info(' Modified Shepard weighting power %5.2f', self.power)
                    self.log.info(' Region of influence radius (pixels) %6.2f', self.roi)

                    result = straylight.correct_mrs_modshepard(input_model,
                                                               regions,
                                                               self.roi,
                                                               self.power)

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
