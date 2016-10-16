#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from .srctype import set_source_type


class SourceTypeStep(Step):
    """
    SourceTypeStep: Selects and sets a source type based on various inputs.
    The source type is used in later calibrations to determine the appropriate
    methods to use. Input comes from either the SRCTYPE keyword value, which
    is populated from user info in the APT, or the NIRSpec MSA planning tool.
    """

    spec = """
    """

    def process(self, input):

        with datamodels.open(input) as input_model:

            # Call the source selection routine
            result = set_source_type(input_model)

            if result is None:
                result = input_model.copy()
                result.meta.cal_step.srctype = 'SKIPPED'
            else:
                result.meta.cal_step.srctype = 'COMPLETE'

        return result


if __name__ == '__main__':
    cmdline.step_script(SourceTypeStep)
