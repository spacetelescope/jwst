import gc
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import reset_sub
from jwst.lib.basic_utils import use_datamodel, copy_datamodel

__all__ = ["ResetStep"]


class ResetStep(Step):
    """
    ResetStep: Performs a reset  correction by subtracting
    the reset correction reference data from the input science data model.
    """

    class_alias = "reset"

    spec = """
    """

    reference_file_types = ['reset']

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model)

        result, input_model = copy_datamodel(input_model, self.parent)

        # check the data is MIRI data
        detector = result.meta.instrument.detector
        if detector.startswith('MIR'):

            # Get the name of the reset reference file to use
            self.reset_name = self.get_reference_file(result, 'reset')
            self.log.info('Using RESET reference file %s', self.reset_name)

            # Check for a valid reference file
            if self.reset_name == 'N/A':
                self.log.warning('No RESET reference file found')
                self.log.warning('Reset step will be skipped')
                result.meta.cal_step.reset = 'SKIPPED'
                return result

            # Open the reset ref file data model
            reset_model = datamodels.ResetModel(self.reset_name)

            # Do the reset correction subtraction
            result = reset_sub.do_correction(result, reset_model)

            # Close the reference file and update the step status
            reset_model.close()
            result.meta.cal_step.reset = 'COMPLETE'

        else:
            self.log.warning('Reset Correction is only for MIRI data')
            self.log.warning('Reset step will be skipped')
            result.meta.cal_step.reset = 'SKIPPED'

        gc.collect()
        return result
