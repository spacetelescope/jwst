#! /usr/bin/env python

from jwst.stpipe import Step
from .. import datamodels
from . import dq_initialization

class DQInitStep(Step):
    """

    DQInitStep:  Initialize the Data Quality extension from the
    mask reference file.  Also initialize the error extension

    """

    reference_file_types = ['mask']

    def process(self, input):

        with datamodels.open(input) as input_model:

            # Retreive the mask reference file name
            self.mask_filename = self.get_reference_file(input_model, 'mask')
            self.log.info('Using MASK reference file %s', self.mask_filename)

            # Check for a valid reference file
            if self.mask_filename == 'N/A':
                self.log.warning('No MASK reference file found')
                self.log.warning('DQ initialization step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.dq_init = 'SKIPPED'
                return result

            # Load the reference file
            mask_model = datamodels.MaskModel(self.mask_filename)

            # Apply the step
            result = dq_initialization.correct_model(input_model, mask_model)

            # Close the reference file
            mask_model.close()

        return result
