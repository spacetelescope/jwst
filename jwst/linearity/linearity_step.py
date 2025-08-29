from stdatamodels.jwst import datamodels

from jwst.linearity import linearity
from jwst.stpipe import Step

__all__ = ["LinearityStep"]


class LinearityStep(Step):
    """Perform a correction for non-linear detector response, using the polynomial method."""

    class_alias = "linearity"

    spec = """
    """  # noqa: E501

    reference_file_types = ["linearity"]

    def process(self, step_input):
        """
        Read in linearity correction and apply it to science data.

        Parameters
        ----------
        step_input : RampModel
            The input ramp model.

        Returns
        -------
        result : RampModel
            The output ramp model with linearity correction applied.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # Work on a copy
            result = input_model.copy()

            # Get the name of the linearity reference file to use
            self.lin_name = self.get_reference_file(result, "linearity")
            self.log.info("Using Linearity reference file %s", self.lin_name)

            # Check for a valid reference file
            if self.lin_name == "N/A":
                self.log.warning("No Linearity reference file found")
                self.log.warning("Linearity step will be skipped")
                result.meta.cal_step.linearity = "SKIPPED"
                return result

            # Open the linearity reference file data model
            lin_model = datamodels.LinearityModel(self.lin_name)

            # Do the linearity correction
            result = linearity.do_correction(result, lin_model)
            result.meta.cal_step.linearity = "COMPLETE"

            # Cleanup
            del lin_model

        return result
