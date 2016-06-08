from jwst.stpipe import Step
from jwst import datamodels
from . import lastframe_sub


class LastFrameStep( Step ):
    """
    LastFrameStep: This is a MIRI specific task to correct the last frame by
    subtracting the array contained in the last frame reference file from the
    last frame.
    """

    reference_file_types = ['lastframe']

    def process(self, input):

        # Open the input data model
        with models.open(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if detector[:3] == 'MIR':

                # Get the name of the reset reference file to use
                self.filename = self.get_reference_file(input_model, 'lastframe')
                self.log.info('Using Last Frame reference file %s', self.filename)

                # Check for a valid reference file
                if self.filename == 'N/A':
                    self.log.warning('No Last Frame reference file found')
                    self.log.warning('Last frame step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.lastframe = 'SKIPPED'
                    return result

                # Open the last frame ref file data model
                lastframe_model = models.LastFrameModel(self.filename)

                # Do the lastframe correction subtraction
                result = lastframe_sub.do_correction(input_model, lastframe_model)

                # Close the reference file and update the step status
                lastframe_model.close()
                result.meta.cal_step.lastframe = 'COMPLETE'

            else:
                self.log.warning('Last Frame Correction is only for MIRI data')
                self.log.warning('Last frame step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.lastframe = 'SKIPPED'

        return result

