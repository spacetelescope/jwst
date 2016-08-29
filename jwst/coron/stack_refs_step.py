#! /usr/bin/env python

from ..stpipe import Step, cmdline

from .. import datamodels
from . import stack_refs

class StackRefsStep(Step):

    """
    StackRefsStep: Stack multiple PSF reference exposures into a
    single CubeModel, for use by subsequent coronagraphic steps.
    """

    spec = """
    """

    def process(self, input):

        # Open the inputs
        with datamodels.open(input) as input_models:

            # Call the stacking routine
            output_model = stack_refs.make_cube(input_models)

        return output_model


if __name__ == '__main__':
    cmdline.step_script(StackRefsStep)
