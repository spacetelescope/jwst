#! /usr/bin/env python
import logging

from jwst.datamodels import ModelLibrary
from jwst.stpipe.utilities import record_step_status

from jwst.stpipe import Step
from .moving_target_wcs import assign_moving_target_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

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
        input_lib : `~jwst.datamodels.ModelLibrary`
            A collection of data models.

        Returns
        -------
        `~jwst.datamodels.ModelLibrary`
            The modified data models.
        """
        if not isinstance(input_lib, ModelLibrary):
            try:
                input_lib = ModelLibrary(input_lib)
            except Exception:
                log.warning("Input data type is not supported.")
                record_step_status(input_lib, "assign_mtwcs", False)
                return input_lib

        result = assign_moving_target_wcs(input_lib)
        return result
