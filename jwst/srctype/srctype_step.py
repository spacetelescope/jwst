#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from .srctype import set_source_type

__all__ = ["SourceTypeStep"]


class SourceTypeStep(Step):
    """
    SourceTypeStep: Selects and sets a source type based on various inputs.
    The source type is used in later calibrations to determine the appropriate
    methods to use. Input comes from either the SRCTYAPT keyword value, which
    is populated from user info in the APT, or the NIRSpec MSA planning tool.
    """

    spec = """
    """

    def process(self, input):

        input_model = datamodels.open(input)

        # Call the source selection routine
        result = set_source_type(input_model)

        # Set the step status in the output model
        if result is None:
            result = input_model
            result.meta.cal_step.srctype = 'SKIPPED'
        else:
            result.meta.cal_step.srctype = 'COMPLETE'

        return result
