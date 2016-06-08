from jwst.stpipe import Step
from jwst import datamodels
from jwst.background import subtract_images


class ImprintStep(Step):
    """
    ImprintStep: Removes NIRSpec MSA imprint structure from an exposure
    by subtracting an imprint exposure.
    """

    spec = """
    """

    def process(self, input, imprint):

        # Open the input data model
        with models.open(input) as input_model:

            # Open the imprint exposure data model
            imprint_model = models.open(imprint)

            # Do the imprint exposure subtraction
            result = subtract_images.subtract(input_model, imprint_model)

            # Update the step status and close the imprint model
            result.meta.cal_step.imprint = 'COMPLETE'
            imprint_model.close()

        return result

