#! /usr/bin/env python

from .. import datamodels
from ..stpipe import Step
from . import combine1d


__all__ = ["Combine1dStep"]


class Combine1dStep(Step):
    """

    Combine1dStep: Combine 1-D spectra

    """

    spec = """
    # integration_time or exposure_time.
    exptime_key = string(default="exposure_time")
    # Interpolation between pixels.  This is currently not used.
    interpolation = string(default="nearest")
    # Set to True if the input spectra are just background.
    background = boolean(default=False)
    """

    def process(self, input_file):

        with datamodels.open(input_file) as input_model:
            result = combine1d.do_combine(input_model,
                                          self.exptime_key,
                                          self.background)

        return result
