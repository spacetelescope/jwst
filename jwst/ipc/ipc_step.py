from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import ipc_corr

__all__ = ["IPCStep"]


class IPCStep(Step):
    """Perform IPC correction by convolving a science datamodel with IPC reference data."""

    class_alias = "ipc"

    spec = """
    """  # noqa: E501

    reference_file_types = ["ipc"]

    def process(self, step_input):
        """
        Apply the IPC correction.

        Parameters
        ----------
        step_input : data model object
            Science data model to be corrected.

        Returns
        -------
        data model object
            IPC-corrected science data model.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # Get the name of the ipc reference file to use
            self.ipc_name = self.get_reference_file(input_model, "ipc")
            self.log.info("Using IPC reference file %s", self.ipc_name)

            # Check for a valid reference file
            if self.ipc_name == "N/A":
                self.log.warning("No IPC reference file found")
                self.log.warning("IPC step will be skipped")
                input_model.meta.cal_step.ipc = "SKIPPED"
                return input_model

            # Open the ipc reference file data model
            ipc_model = datamodels.IPCModel(self.ipc_name)

            # Work on a copy
            result = input_model.copy()

            # Do the ipc correction
            result = ipc_corr.do_correction(result, ipc_model)
            result.meta.cal_step.ipc = "COMPLETE"

            # Cleanup
            del ipc_model

        return result
