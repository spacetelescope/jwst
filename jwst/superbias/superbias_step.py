from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import bias_sub

__all__ = ["SuperBiasStep"]


class SuperBiasStep(Step):
    """Subtract super-bias reference data from the input science data model."""

    class_alias = "superbias"

    spec = """
    """  # noqa: E501

    reference_file_types = ["superbias"]

    def process(self, step_input):
        """
        Apply the superbias correction.

        Parameters
        ----------
        step_input : str or stdatamodels.jwst.datamodels.ramp.RampModel
            Either the path to the file or the science data model to be corrected.

        Returns
        -------
        stdatamodels.jwst.datamodels.ramp.RampModel
            Superbias-corrected science data model.
        """
        # Open the input data model
        with datamodels.open(step_input) as input_model:
            # Get the name of the superbias reference file to use
            self.bias_name = self.get_reference_file(input_model, "superbias")
            self.log.info("Using SUPERBIAS reference file %s", self.bias_name)

            # Check for a valid reference file
            if self.bias_name == "N/A":
                self.log.warning("No SUPERBIAS reference file found")
                self.log.warning("Superbias step will be skipped")
                input_model.meta.cal_step.superbias = "SKIPPED"
                return input_model

            # Open the superbias ref file data model
            bias_model = datamodels.SuperBiasModel(self.bias_name)

            # Work on a copy
            result = input_model.copy()

            # Do the bias subtraction
            result = bias_sub.do_correction(result, bias_model)
            result.meta.cal_step.superbias = "COMPLETE"

            # Cleanup
            del bias_model

        return result
