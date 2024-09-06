#! /usr/bin/env python
import logging

from jwst.datamodels import ModelLibrary
from jwst.stpipe.utilities import record_step_status

from ..stpipe import Step
from .moving_target_wcs import assign_moving_target_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["AssignMTWcsStep"]


class AssignMTWcsStep(Step):
    """
    AssignMTWcsStep: Create a gWCS object for a moving target.

    Parameters
    ----------
    input : `~jwst.associations.Association`
        A JWST association file.
    """

    class_alias = "assign_mtwcs"

    spec = """
        suffix = string(default='assign_mtwcs')    # Default suffix for output files
        output_use_model = boolean(default=True)   # When saving use `DataModel.meta.filename`
    """

    def process(self, input):

        if not isinstance(input, ModelLibrary):
            try:
                input = ModelLibrary(input)
            except Exception:
                log.warning("Input data type is not supported.")
                record_step_status(input, "assign_mtwcs", False)
                return input
            
        result = assign_moving_target_wcs(input)
        return result
