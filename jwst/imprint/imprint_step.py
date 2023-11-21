from stdatamodels.jwst import datamodels

from ..stpipe import Step
from ..background import subtract_images

__all__ = ["ImprintStep"]


class ImprintStep(Step):
    """
    ImprintStep: Removes NIRSpec MSA imprint structure from an exposure
    by subtracting an imprint (a.k.a. leakcal) exposure.
    """

    class_alias = "imprint"

    spec = """
    """

    def process(self, input, imprint):

        # Open the input science image and get its dither pattern position number
        input_model = datamodels.open(input)
        pos_no = input_model.meta.dither.position_number
        obs_no = input_model.meta.observation.observation_number

        # If there is only one imprint image listed in the association,
        # use it for all science images.
        # If there is more than one imprint image in the association,
        # then select the imprint image to be used by matching the dither pattern position
        # number and observation number to the science image.
        num_imprint = len(imprint)
        if num_imprint == 1:
            match = 0  # use the first, and only, imprint image
        else:
            match = None
            for i in range(num_imprint):
                # open imprint image and read in pos_no to find match
                imprint_model = datamodels.open(imprint[i])
                imprint_pos_no = imprint_model.meta.dither.position_number
                imprint_obs_no = imprint_model.meta.observation.observation_number
                imprint_model.close()
                if pos_no == imprint_pos_no and obs_no == imprint_obs_no:
                    match = i
                    break

        # Copy the input image to the output (just in case no matching imprint image was found)
        result = input_model.copy()

        if match is not None:
            # Subtract the matching imprint image
            imprint_model = datamodels.open(imprint[match])
            self.log.info(f"Subtracting imprint image {imprint_model.meta.filename}")
            result = subtract_images.subtract(input_model, imprint_model)

            # Update the step status and close the imprint model
            result.meta.cal_step.imprint = 'COMPLETE'
            imprint_model.close()
        else:
            self.log.warning(f'No matching imprint image was found for {input}')
            self.log.warning('Step will be skipped')
            result.meta.cal_step.imprint = 'SKIPPED'

        input_model.close()
        return result
