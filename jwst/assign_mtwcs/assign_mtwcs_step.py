#! /usr/bin/env python
import logging
from ..stpipe import Step
from .. import datamodels
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
        if isinstance(input, str):
            input = datamodels.open(input)

        # Can't apply the step if we aren't given a ModelContainer as input
        if not isinstance(input, datamodels.ModelContainer):
            log.warning("Input data type is not supported.")
            # raise ValueError("Expected input to be an association file name or a ModelContainer.")
            input.meta.cal_step.assign_mtwcs = 'SKIPPED'
            return input

        # Apply the step
        result = assign_moving_target_wcs(input)

        return result
