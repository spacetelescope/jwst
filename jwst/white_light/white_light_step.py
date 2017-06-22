#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from .white_light import white_light


class WhiteLightStep(Step):
    """
    WhiteLightStep: Computes integrated flux as a function of time for a
    multi-integration spectroscopic observation.
    """

    spec = """
    """

    def process(self, input):

        # Load the input
        with datamodels.open(input) as input_model:

            # Call the white light curve generation routine
            result = white_light(input_model)

            # Write the output catalog
            index = input_model.meta.filename.rfind('_')
            cat_name = input_model.meta.filename[:index + 1] + 'cat.ecsv'
            result.write(cat_name, format='ascii.ecsv')

        return
