#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from .white_light import white_light

__all__ = ["WhiteLightStep"]


class WhiteLightStep(Step):
    """
    WhiteLightStep: Computes integrated flux as a function of time for a
    multi-integration spectroscopic observation.
    """

    spec = """
    min_wavelength     = float(default=None)      # Default wavelength minimum for integration
    max_wavelength     = float(default=None)      # Default wavelength maximum for integration
    output_ext         = string(default='.ecsv')  # Default type of output
    suffix             = string(default='whtlt')  # Default suffix for output files
    """

    def process(self, input):

        # Load the input
        with datamodels.open(input) as input_model:

            # Call the white light curve generation routine
            result = white_light(input_model, self.min_wavelength, self.max_wavelength)

            # Write the output catalog
            if self.save_results:
                output_path = self.make_output_path()
                result.write(output_path, format='ascii.ecsv', overwrite=True)

        return result
