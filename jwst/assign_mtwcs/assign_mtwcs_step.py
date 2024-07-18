#! /usr/bin/env python
import logging

from jwst.datamodels import ModelLibrary, ModelContainer
from jwst.datamodels.library import container_to_library, library_to_container

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
        if isinstance(input, (str, dict)):
            input = ModelLibrary(input)
        elif isinstance(input, ModelContainer):
            input = container_to_library(input)

        # Can't apply the step if we aren't given a ModelLibrary as input
        if not isinstance(input, ModelLibrary):
            log.warning("Input data type is not supported.")
            # raise ValueError("Expected input to be an association file name or a ModelContainer.")
            input.meta.cal_step.assign_mtwcs = 'SKIPPED'
            return library_to_container(input)

        # Apply the step
        result = assign_moving_target_wcs(input)

        return library_to_container(result)
