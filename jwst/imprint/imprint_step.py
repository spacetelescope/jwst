from stdatamodels.jwst import datamodels

from ..stpipe import Step
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

    def process(self, input, imprint):

        # subtract leakcal (imprint) image
        # If there is only one imprint image is in the association use it  for all the data.
        # If there is more than one imprint image in the association then select
        # the imprint  image to subtract based on position number (DataModel.meta.dither.position_number) &
        # observation number (Datamodel.meta.observation.observation_number)

        # Open the input data model and get position number of image
        input_model = datamodels.open(input)
        pos_no = input_model.meta.dither.position_number
        obs_no = input_model.meta.observation.observation_number

        # find imprint that goes with input image - if there is only 1 imprint image - use it.
        num_imprint = len(imprint)
        match = None
        if num_imprint == 1:
            match = 0
        else:
            i = 0
            while match is None and i < num_imprint:
                # open imprint images and read in pos_no to find match
                imprint_model = datamodels.open(imprint[i])
                imprint_pos_no = imprint_model.meta.dither.position_number
                imprint_obs_no = imprint_model.meta.observation.observation_number
                imprint_model.close()
                if pos_no == imprint_pos_no and obs_no == imprint_obs_no:
                    match = i
                i = i + 1

        # Initialize result to be input (just in case no matching imprint image was found)
        result = input_model.copy()

        if match is not None:
            imprint_model = datamodels.open(imprint[match])
            result = subtract_images.subtract(input_model, imprint_model)

            # Update the step status and close the imprint model
            result.meta.cal_step.imprint = 'COMPLETE'
            imprint_model.close()
        else:
            self.log.info(f'No imprint image was found for {input}')
            result.meta.cal_step.imprint = 'SKIPPED'

        input_model.close()
        return result
