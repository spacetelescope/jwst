from ..stpipe import Step
from .. import datamodels
from . import gain_scale


__all__ = ["GainScaleStep"]


class GainScaleStep(Step):
    """
    GainScaleStep: Rescales countrate data to account for use of a
    non-standard gain value. All integrations are multiplied by the
    factor GAINFACT.
    """

    class_alias = "gain_scale"

    reference_file_types = ['gain']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # Is the gain_factor already populated in the input model?
            if input_model.meta.exposure.gain_factor is None:

                # Try to get the gain factor from the gain reference file
                self.gain_filename = self.get_reference_file(input_model, 'gain')
                gain_model = datamodels.GainModel(self.gain_filename)

                # Try to read the GAINFACT keyword value
                if gain_model.meta.exposure.gain_factor is None:
                    self.log.info('GAINFACT not found in gain reference file')
                    self.log.info('Step will be skipped')
                    input_model.meta.cal_step.gain_scale = 'SKIPPED'
                    gain_model.close()
                    return input_model
                else:
                    gain_factor = gain_model.meta.exposure.gain_factor
                    gain_model.close()

            else:
                gain_factor = input_model.meta.exposure.gain_factor

            # Do the scaling
            result = gain_scale.do_correction(input_model, gain_factor)

        return result
