#! /usr/bin/env python

import os
from astropy.io import fits             # xxx temporary
from ..stpipe import Step
from .. import datamodels
from . import extract

class Extract1dStep(Step):
    """
    Extract1dStep: Extract a 1-d spectrum from 2-d data
    """

    spec = """
    # Boxcar smoothing width for background regions.
    smoothing_length = integer(default=None)
    # Order of polynomial fit to one column (or row if the dispersion
    # direction is vertical) of background regions.
    bkg_order = integer(default=None, min=0)
    """

    reference_file_types = ['extract1d']

    def process(self, input):

        # Open the input and figure out what type of model it is
        input_model = datamodels.open(input)
        data_model_from_header = input_model.meta.model_type
        self.log.debug("Data model from header = %s", data_model_from_header)

        if isinstance(input_model, datamodels.CubeModel):
            # It's a 3-D multi-integration model
            self.log.debug('Input is a CubeModel for a multiple integ. file')
        elif isinstance(input_model, datamodels.ImageModel):
            # It's a single 2-D image
            self.log.debug('Input is an ImageModel')
        elif isinstance(input_model, datamodels.MultiSlitModel):
            self.log.debug('Input is a MultiSlitModel')
        elif isinstance(input_model, datamodels.IFUCubeModel):
            self.log.debug('Input is an IFUCubeModel')
        elif isinstance(input_model, datamodels.DrizProductModel):
            # Resampled 2-D data
            self.log.debug('Input is a DrizProductModel')
        else:
            self.log.warning('Input is a %s,', str(type(input_model)))
            self.log.warning('which was not expected for extract_1d.')

        # Get the reference file name
        self.ref_file = self.get_reference_file(input_model, 'extract1d')
        self.log.info('Using EXTRACT1D reference file %s', self.ref_file)

        # Do the extraction
        result = extract.do_extract1d(input_model, self.ref_file,
                                      self.smoothing_length, self.bkg_order)

        # Set the step flag to complete
        result.meta.cal_step.extract_1d = 'COMPLETE'

        return result

def extract_1d_correction(input):
    a = Extract1dStep()
    result = a.process(input)
    return result
