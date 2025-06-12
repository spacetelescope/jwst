from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from . import gain_scale

__all__ = ["GainScaleStep"]


class GainScaleStep(Step):
    """
    Rescale all integrations in an exposure by gain_factor.

    GainScaleStep: Rescales countrate data to account for use of a
    non-standard gain value. All integrations are multiplied by the
    factor GAINFACT.
    """

    class_alias = "gain_scale"

    spec = """
    """  # noqa: E501
    reference_file_types = ["gain"]

    def process(self, step_input):
        """
        Perform gain scale step.

        Parameters
        ----------
        step_input : datamodel
            Input datamodel on which to perform gain scale step.

        Returns
        -------
        result : datamodel
            Output datamodel on which the gain scale step has been performed.
        """
        # Open the input data model
        with datamodels.open(step_input) as input_model:
            # Is the gain_factor already populated in the input model?
            if input_model.meta.exposure.gain_factor is None:
                # Try to get the gain factor from the gain reference file
                self.gain_filename = self.get_reference_file(input_model, "gain")
                gain_model = datamodels.GainModel(self.gain_filename)

                # Try to read the GAINFACT keyword value
                if gain_model.meta.exposure.gain_factor is None:
                    self.log.info("GAINFACT not found in gain reference file")
                    self.log.info("Step will be skipped")
                    input_model.meta.cal_step.gain_scale = "SKIPPED"
                    del gain_model
                    return input_model
                else:
                    gain_factor = gain_model.meta.exposure.gain_factor
                    del gain_model

            else:
                gain_factor = input_model.meta.exposure.gain_factor

            # Work on a copy
            result = input_model.copy()

            # Do the scaling
            result = gain_scale.do_correction(result, gain_factor)

        return result
