import logging

from stdatamodels.jwst import datamodels

from jwst.ipc import ipc_corr
from jwst.stpipe import Step

__all__ = ["IPCStep"]

log = logging.getLogger(__name__)


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
        step_input : `~stdatamodels.jwst.datamodels.RampModel`
            Science data model to be corrected.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            IPC-corrected science data model.
        """
        # Open the input data model
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        # Get the name of the ipc reference file to use
        ipc_name = self.get_reference_file(result, "ipc")
        log.info("Using IPC reference file %s", ipc_name)

        # Check for a valid reference file
        if ipc_name == "N/A":
            log.warning("No IPC reference file found")
            log.warning("IPC step will be skipped")
            result.meta.cal_step.ipc = "SKIPPED"
            return result

        # Open the ipc reference file data model
        ipc_model = datamodels.IPCModel(ipc_name)

        # Do the ipc correction
        result = ipc_corr.do_correction(result, ipc_model)
        result.meta.cal_step.ipc = "COMPLETE"

        # Cleanup
        del ipc_model

        return result
