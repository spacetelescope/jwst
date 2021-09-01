from ..stpipe import Step
from .. import datamodels
from . import reset_sub

__all__ = ["ResetStep"]


class ResetStep(Step):
    """
    ResetStep: Performs a reset  correction by subtracting
    the reset correction reference data from the input science data model.
    """

    class_alias = "reset"

    reference_file_types = ['reset']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if detector.startswith('MIR'):

                # Get the name of the reset reference file to use
                self.reset_name = self.get_reference_file(input_model, 'reset')
                self.log.info('Using RESET reference file %s', self.reset_name)

                # Check for a valid reference file
                if self.reset_name == 'N/A':
                    self.log.warning('No RESET reference file found')
                    self.log.warning('Reset step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.reset = 'SKIPPED'
                    return result

                # Open the reset ref file data model
                reset_model = datamodels.ResetModel(self.reset_name)

                # Do the reset correction subtraction
                result = reset_sub.do_correction(input_model, reset_model)

                # Close the reference file and update the step status
                reset_model.close()
                result.meta.cal_step.reset = 'COMPLETE'

            else:
                self.log.warning('Reset Correction is only for MIRI data')
                self.log.warning('Reset step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.reset = 'SKIPPED'

        return result
