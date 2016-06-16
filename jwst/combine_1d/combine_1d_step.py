#! /usr/bin/env python

import time
import os
from ..stpipe import Step
from . import combine1d

class Combine1dStep(Step):
    """

    Combine1dStep: Combine 1-D spectra

    """

    spec = """
    # integration_time or exposure_time.
    exptime_key = string(default="integration_time")
    # Interpolation between pixels.
    interpolation = string(default="nearest")
    """

    def process(self, input_file):

        combine1d.correct_model(input_file,
                                self.exptime_key,
                                self.interpolation)

def combine_1d_correction(input):
    a = Combine1dStep()
    a.process(input)
