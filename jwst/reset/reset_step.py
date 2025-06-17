from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import reset_sub

__all__ = ["ResetStep"]


class ResetStep(Step):
    """Subtract the reset correction reference data from the MIRI input science data model."""

    class_alias = "reset"

    spec = """
    """  # noqa: E501

    reference_file_types = ["reset"]

    def process(self, step_input):
        """
        Subtract the reset correction from the MIRI ramp model.

        Parameters
        ----------
        step_input : DataModel
           Input datamodel to be corrected

        Returns
        -------
        reset: DataModel
           The reset corrected ramp model.
        """
        # Open the input data model
        with datamodels.open(step_input) as input_model:
            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if not detector.startswith("MIR"):
                self.log.warning("Reset Correction is only for MIRI data")
                self.log.warning("Reset step will be skipped")
                input_model.meta.cal_step.reset = "SKIPPED"
                return input_model

            # Get the name of the reset reference file to use
            self.reset_name = self.get_reference_file(input_model, "reset")
            self.log.info(f"Using RESET reference file {self.reset_name}")

            # Check for a valid reference file
            if self.reset_name == "N/A":
                self.log.warning("No RESET reference file found")
                self.log.warning("Reset step will be skipped")
                input_model.meta.cal_step.reset = "SKIPPED"
                return input_model

            # Open the reset ref file data model
            reset_model = datamodels.ResetModel(self.reset_name)

            # Work on a copy
            result = input_model.copy()

            # Do the reset correction subtraction
            result = reset_sub.do_correction(result, reset_model)
            result.meta.cal_step.reset = "COMPLETE"

            # Cleanup
            del reset_model

        return result
