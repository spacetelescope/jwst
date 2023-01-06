from ..stpipe import Step
from .. import datamodels
from ..background import subtract_images

__all__ = ["ImprintStep"]


class ImprintStep(Step):
    """
    ImprintStep: Removes NIRSpec MSA imprint structure from an exposure
    by subtracting an imprint exposure.
    """

    class_alias = "imprint"

    spec = """
    """

    def process(self, input, imprint, imprint_pos_no):

        # imprint_pos_no - a list of the position numbers corresponding to the list of imprints
        
        # subtract leakcal (imprint) image
        # If only 1 imprint image is in the association use for for science and if there a background
        # If more than 1 imprint image exists in the association then select the imprint  image to
        # subtract based on position number.

        # Open the input data model and get position number of image
        input_model = datamodels.open(input)
        pos_no = input_model.meta.dither.position_number

        # find imprint that goes with inpyt image - if there is only 1 imprint image - use it.
        num_imprint = len(imprint)
        match = None
        if num_imprint == 1:
            match = 0
        else:
            i = 0
            while match is None and i < num_imprint:
                if pos_no == imprint_pos_no[i]:
                    match = i
                i = i + 1

        # Initialize result to be input (just in case not matching imprint image was found)
        result = input_model.copy()

        if match is not None:
            imprint_model = datamodels.open(imprint[match])
            result = subtract_images.subtract(input_model, imprint_model)

            # Update the step status and close the imprint model
            result.meta.cal_step.imprint = 'COMPLETE'
            input_model.close()
            imprint_model.close()
        else:
            log.info(f'No imprint image was found for {input}')
        return result
