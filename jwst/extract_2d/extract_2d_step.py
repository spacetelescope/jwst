#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import extract_2d

class Extract2dStep(Step):
    """
    This Step performs a 2D extraction of spectra.
    """

    spec = """
        which_subarray = string(default = None)
    """

    def process(self, input_file):

        with datamodels.open(input_file) as dm:

            output_model = extract_2d.extract2d(dm, self.which_subarray)
        
        return output_model


if __name__ == '__main__':
    cmdline.step_script(extract_2d_step)
