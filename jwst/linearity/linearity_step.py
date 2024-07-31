import gc
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import linearity
from jwst.lib.basic_utils import use_datamodel, copy_datamodel

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

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model, model_class=datamodels.RampModel)

        result, input_model = copy_datamodel(input_model, self.parent)

        # Get the name of the linearity reference file to use
        self.lin_name = self.get_reference_file(result, 'linearity')
        self.log.info('Using Linearity reference file %s', self.lin_name)

        # Check for a valid reference file
        if self.lin_name == 'N/A':
            self.log.warning('No Linearity reference file found')
            self.log.warning('Linearity step will be skipped')
            result.meta.cal_step.linearity = 'SKIPPED'
            gc.collect()
            return result

        # Open the linearity reference file data model
        lin_model = datamodels.LinearityModel(self.lin_name)

        # Do the linearity correction
        result = linearity.do_correction(result, lin_model)

        # Close the reference file and update the step status
        del lin_model
        result.meta.cal_step.linearity = 'COMPLETE'

        gc.collect()
        return result
