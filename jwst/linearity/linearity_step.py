from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import linearity

__all__ = ["LinearityStep"]


class LinearityStep(Step):
    """
    LinearityStep: This step performs a correction for non-linear
    detector response, using the "classic" polynomial method.
    """

    class_alias = "linearity"

    spec = """
    """

    reference_file_types = ['linearity']

    def process(self, step_input):

        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:

            # Get the name of the linearity reference file to use
            self.lin_name = self.get_reference_file(input_model, 'linearity')
            self.log.info('Using Linearity reference file %s', self.lin_name)

            # Check for a valid reference file
            if self.lin_name == 'N/A':
                self.log.warning('No Linearity reference file found')
                self.log.warning('Linearity step will be skipped')
                input_model.meta.cal_step.linearity = 'SKIPPED'
                return input_model

            # Open the linearity reference file data model
            lin_model = datamodels.LinearityModel(self.lin_name)

            # Work on a copy
            result = input_model.copy()

            # Do the linearity correction
            result = linearity.do_correction(result, lin_model)
            result.meta.cal_step.linearity = 'COMPLETE'

            # Cleanup
            del lin_model

        return result
