from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import nsclean

__all__ = ["NSCleanStep"]


class NSCleanStep(Step):
    """
    NSCleanStep: This step performs 1/f noise correction ("cleaning")
    of NIRSpec images, using the "NSClean" method.
    """

    class_alias = "nsclean"

    spec = """
        n_sigma = float(default=5.0)
    """

    #reference_file_types = ['nsc_mask']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # Do the NSClean correction
            result = nsclean.do_correction(input_model, self.n_sigma)

            # Update the step status
            result.meta.cal_step.nsclean = 'COMPLETE'

        return result
