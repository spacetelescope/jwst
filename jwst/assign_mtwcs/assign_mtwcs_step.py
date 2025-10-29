#! /usr/bin/env python
import logging

from jwst.assign_mtwcs.moving_target_wcs import assign_moving_target_wcs
from jwst.datamodels import ModelLibrary
from jwst.stpipe import Step
from jwst.stpipe.utilities import record_step_status

log = logging.getLogger(__name__)

__all__ = ["AssignMTWcsStep"]


class AssignMTWcsStep(Step):
    """Create a gWCS object for a moving target."""

    class_alias = "assign_mtwcs"

    spec = """
        suffix = string(default='assign_mtwcs')    # Default suffix for output files
        output_use_model = boolean(default=True)   # When saving use `DataModel.meta.filename`
    """  # noqa: E501

    def process(self, input_lib):
        """
        Run the assign_mtwcs step.

        Parameters
        ----------
        input_lib : `~jwst.datamodels.library.ModelLibrary`
            A collection of data models.

        Returns
        -------
        `~jwst.datamodels.library.ModelLibrary`
            The modified data models.
        """
        # Open the input, making a copy if necessary
        output_lib = self.prepare_output(input_lib)
        if not isinstance(output_lib, ModelLibrary):
            try:
                output_lib = ModelLibrary(output_lib)
            except (ValueError, TypeError) as err:
                log.warning("Input data type is not supported.")
                log.debug(f"Error was: {err}")
                record_step_status(output_lib, "assign_mtwcs", False)
                return output_lib

        result = assign_moving_target_wcs(output_lib)
        record_step_status(result, "assign_mtwcs", True)
        return result
