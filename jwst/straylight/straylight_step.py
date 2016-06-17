#! /usr/bin/env python

from jwst.stpipe import Step, cmdline
from jwst import datamodels
from . import straylight


class StraylightStep (Step):
    """
    StraylightStep: Performs straylight correction image using a Mask file.
    """
    reference_file_types = ['straymask']

    def process(self, input):

        # Open the input data model
        with models.ImageModel(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if detector == 'MIRIFUSHORT':

                # Get the name of the straylight reference file
                self.straylight_name = self.get_reference_file(input_model,
                                                               'straymask')
                self.log.info('Using STRAYLIGHT reference file %s',
                                                      self.straylight_name)

                # Check for a valid reference file
                if self.straylight_name == 'N/A':
                    self.log.warning('No STRAYLIGHT reference file found')
                    self.log.warning('Straylight step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.straylight = 'SKIPPED'
                    return result

                # Open the straylight mask ref file data model
                straylight_model = models.StrayLightModel(self.straylight_name)

                # Do the correction
                result = straylight.correct_MRS(input_model, straylight_model)

                # Close the reference file and update the step status
                straylight_model.close()
                result.meta.cal_step.straylight = 'COMPLETE'

            else:
                self.log.warning('Straylight correction not defined for detector %s',
                                  detector)
                self.log.warning('Straylight step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.straylight = 'SKIPPED'

        return result


if __name__ == '__main__':
    cmdline.step_script(straylight_step)
