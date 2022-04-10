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

            # Check for a valid reference file
            # Use the Regions reference file set to 20% throughput threshold
            #self.straylight_name = self.get_reference_file(input_model,
            #                                               'regions')
            # Read in the hacked new reference file
            # (to be replace with a proper ref file read from CRDS)
            self.mrsxart_name = os.path.join(os.path.expandvars('$MIRI3D_DATA_DIR'),'mrsxartcor/temp/miri-mrsxartcor-59663.fits')

            self.log.info('Using mrsxart reference file %s',
                          self.mrsxart_name)

            if self.mrsxart_name == 'N/A':
                self.log.warning('No MRSXART reference file found')
                self.log.warning('Straylight step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.straylight = 'SKIPPED'
                return result

            # Read in the reference file
            #allregions = datamodels.RegionsModel(self.straylight_name)
            # Use 20% throughput array
            #regions = (allregions.regions)[2, :, :].copy()
            #self.log.info(' Using 20% throughput threshold.')
            #allregions.close()
            modelpars = fits.open(self.mrsxart_name)

            # Apply the correction
            result = straylight.correct_xartifact(input_model, modelpars)

            result.meta.cal_step.straylight = 'COMPLETE'

        return result


class ErrorNoAssignWCS(Exception):
    pass
