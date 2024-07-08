import gc
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import ipc_corr
from jwst.lib.basic_utils import use_datamodel, copy_datamodel

__all__ = ["IPCStep"]


class IPCStep(Step):
    """
    IPCStep: Performs IPC correction by convolving the input science
    data model with the IPC reference data.
    """

    class_alias = "ipc"

    spec = """
    """

    reference_file_types = ['ipc']

    def process(self, input_model):
        """Apply the IPC correction.

        Parameters
        ----------
        input : data model object
            Science data model to be corrected.

        Returns
        -------
        data model object
            IPC-corrected science data model.
        """

        # Open the input data model
        input_model = use_datamodel(input_model, model_class=datamodels.RampModel)

        result, input_model = copy_datamodel(input_model, self.parent)

        # Get the name of the ipc reference file to use
        self.ipc_name = self.get_reference_file(result, 'ipc')
        self.log.info('Using IPC reference file %s', self.ipc_name)

        # Check for a valid reference file
        if self.ipc_name == 'N/A':
            self.log.warning('No IPC reference file found')
            self.log.warning('IPC step will be skipped')
            result.meta.cal_step.ipc = 'SKIPPED'
            return result

        # Open the ipc reference file data model
        ipc_model = datamodels.IPCModel(self.ipc_name)

        # Do the ipc correction
        result = ipc_corr.do_correction(result, ipc_model)

        # Close the reference file and update the step status
        del ipc_model
        result.meta.cal_step.ipc = 'COMPLETE'

        gc.collect()
        return result
