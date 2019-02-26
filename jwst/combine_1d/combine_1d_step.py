#! /usr/bin/env python

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
    # Interpolation between pixels.
    interpolation = string(default="nearest")
    # Set to True if the input spectra are just background.
    background = boolean(default=False)
    """

    def process(self, input_file):

        result = combine1d.do_correction(input_file,
                                         self.exptime_key,
                                         self.interpolation,
                                         self.background)

        return result
