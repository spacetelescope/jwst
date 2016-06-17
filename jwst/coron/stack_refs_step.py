#! /usr/bin/env python

from ..stpipe import Step, cmdline

import json
from . import stack_refs

class StackRefsStep(Step):

    """
    StackRefsStep: Stack multiple PSF reference images into a CubeModel, for 
    use by subsequent coronagraphic processing tasks.
    """

    spec = """
    """

    def process(self, input):

        # Load the input list of PSF reference images
        input_table = LoadPsfRefs(input)

        # Call the stacking routine
        cube_model = stack_refs.make_cube(input_table)

        #result.meta.cal_step.psfcube = 'COMPLETE'

        return cube_model


class LoadPsfRefs(object):
    """Class to handle reading a coronagraphic asn table and loading the
       input and output member info into a table model.
    """

    def __init__(self, input):
        self.input = input # keep a record of original input name for later

        if isinstance(input, str):
            self.asn_table = json.load(open(input, 'r'))
        else:
            raise TypeError

        # Extract some values for ease of use
        self.input_filenames = self.get_inputs()
        self.output_filename = self.get_output()

    def get_inputs(self, product=0):
        members = []
        for p in self.asn_table['products'][product]['members']:
            members.append(p['expname'])
        return members

    def get_output(self, product=0):
        return self.asn_table['products'][product]['name']


if __name__ == '__main__':
    cmdline.step_script(StackRefsStep)
