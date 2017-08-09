from ..stpipe import Step
from .. import datamodels
from . import gain_scale


class GainScaleStep(Step):
    """
    GainScaleStep: Rescales countrate data to account for use of a
    non-standard gain value. All integrations are multiplied by the
    factor GAINFACT.
    """

    reference_file_types = ['gain']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # Retrieve the gain reference file name
            self.gain_filename = self.get_reference_file(input_model, 'gain')

            # Load the gain reference file
            gain_model = datamodels.GainModel(self.gain_filename)

            # Try to read the GAINFACT keyword value
            try:
                gain_factor = gain_model.meta.exposure.gain_factor
                gain_model.close()
            except:
                self.log.warning('GAINFACT not found in gain reference file')
                self.log.warning('Step will be skipped')
                input_model.meta.cal_step.gain_scale = 'SKIPPED'
                gain_model.close()
                return input_model

            # Do the scaling
            result = gain_scale.do_correction(input_model, gain_factor)

        return result
