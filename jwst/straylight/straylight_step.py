#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import straylight
from ..datamodels import RegionsModel

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
    reference_file_types = ['straymask']
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
                    # going to use Regions file that is in the ASDF extension
                    assign_wcs = input_model.meta.cal_step.assign_wcs
                    if(assign_wcs != 'COMPLETE'):
                        self.log.warning('Assign_WCS was not run on file, we  need the information of the slice gap locations')
                        raise ErrorNoAssignWCS("Assign WCS has not been run on file")

                    det2ab = input_model.meta.wcs.get_transform('detector', 'alpha_beta')
                    #det2ab is a RegionsSelector model
                    slices = det2ab.label_mapper.mapper

                    self.log.info(' Region of influence radius (pixels) %6.2f', self.roi)
                    self.log.info(' Modified Shepard weighting power %5.2f', self.power)
                # Do the correction
                    result = straylight.correct_MRS_ModShepard(input_model,
                                                               slices,
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

class ErrorNoAssignWCS(Exception):
    pass
