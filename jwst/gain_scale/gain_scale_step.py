from stdatamodels.jwst import datamodels

from jwst.gain_scale import gain_scale
from jwst.stpipe import Step

__all__ = ["GainScaleStep"]


class GainScaleStep(Step):
    """Rescale all integrations in an exposure by a gain factor."""

    class_alias = "gain_scale"

    spec = """
    """  # noqa: E501
    reference_file_types = ["gain"]

    def process(self, step_input):
        """
        Perform gain scale step.

        Rescales countrate data to account for use of a non-standard
        gain value. All integrations are multiplied by the factor GAINFACT.

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
            # Work on a copy
            result = input_model.copy()

            # Is the gain_factor already populated in the input model?
            if result.meta.exposure.gain_factor is None:
                # Try to get the gain factor from the gain reference file
                gain_filename = self.get_reference_file(result, "gain")
                if gain_filename != "N/A":
                    self.log.info("Using GAIN reference file: %s", gain_filename)
                    with datamodels.GainModel(gain_filename) as gain_model:
                        gain_factor = gain_model.meta.exposure.gain_factor
                else:
                    gain_factor = None

                # Try to read the GAINFACT keyword value
                if gain_factor is None:
                    self.log.info("GAINFACT not found in gain reference file")
                    self.log.info("Step will be skipped")
                    result.meta.cal_step.gain_scale = "SKIPPED"
                    return result

            else:
                gain_factor = input_model.meta.exposure.gain_factor

            # Do the scaling
            result = gain_scale.do_correction(result, gain_factor)

        return result
