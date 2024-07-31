from stdatamodels.jwst import datamodels
from ..stpipe import Step
from . import gain_scale
from jwst.lib.basic_utils import use_datamodel, copy_datamodel
import gc

__all__ = ["GainScaleStep"]


class GainScaleStep(Step):
    """
    GainScaleStep: Rescales countrate data to account for use of a
    non-standard gain value. All integrations are multiplied by the
    factor GAINFACT.
    """

    class_alias = "gain_scale"

    spec = """
    """
    reference_file_types = ['gain']

    def process(self, input_model):

        # Open the input data model
        input_model = use_datamodel(input_model)

        result, input_model = copy_datamodel(input_model, self.parent)

        # Is the gain_factor already populated in the input model?
        if result.meta.exposure.gain_factor is None:

            # Try to get the gain factor from the gain reference file
            self.gain_filename = self.get_reference_file(result, 'gain')
            gain_model = datamodels.GainModel(self.gain_filename)

            # Try to read the GAINFACT keyword value
            if gain_model.meta.exposure.gain_factor is None:
                self.log.info('GAINFACT not found in gain reference file')
                self.log.info('Step will be skipped')
                result.meta.cal_step.gain_scale = 'SKIPPED'
                del gain_model
                gc.collect()
                return result
            else:
                gain_factor = gain_model.meta.exposure.gain_factor
                del gain_model

        else:
            gain_factor = result.meta.exposure.gain_factor

        # Do the scaling
        result = gain_scale.do_correction(result, gain_factor)

        gc.collect()
        return result
