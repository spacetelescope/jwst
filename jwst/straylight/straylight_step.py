#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import straylight
from astropy.io import fits
import os

__all__ = ["StraylightStep"]


class StraylightStep (Step):
    """
    StraylightStep: Performs straylight correction image using a Mask file.
    """

    class_alias = "straylight"

    spec = """
         method = option('Nearest','ModShepard',default='ModShepard') #Algorithm method
         roi = integer(2, 1024, default = 50) # Region of interest given as even integer
         power = float(0.1, 5, default = 1.0) # Power of weighting function

    """

    # Change this line
    reference_file_types = ['regions']

    def process(self, input):
        # Open the input data model
        with datamodels.IFUImageModel(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector

            # Check for a valid mrsxartcorr reference file
            self.straylight_name = self.get_reference_file(input_model, 'mrsxartcorr')

            self.log.info('Using mrsxartcorr reference file %s',
                          self.straylight_name)

            if self.straylight_name == 'N/A':
                self.log.warning('No MRSXARTCORR reference file found')
                self.log.warning('Straylight step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.straylight = 'SKIPPED'
                return result

            modelpars = datamodels.MirMrsXArtCorrModel(self.straylight_name)

            # Apply the correction
            result = straylight.correct_xartifact(input_model, modelpars)

            modelpars.close()
            result.meta.cal_step.straylight = 'COMPLETE'

        return result


class ErrorNoAssignWCS(Exception):
    pass
