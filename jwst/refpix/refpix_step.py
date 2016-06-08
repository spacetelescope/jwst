#! /usr/bin/env python

from jwst.stpipe import Step
from jwst import datamodels
from . import reference_pixels
from . import irs2_subtract_reference

class RefPixStep(Step):
    """

    RefPixStep: Use reference pixels to correct bias drifts

    """

    spec = """
        odd_even_columns = boolean(default=True)
        use_side_ref_pixels = boolean(default=True)
        side_smoothing_length = integer(default=11)
        side_gain = float(default=1.0)
        odd_even_rows = boolean(default=True)
    """

    reference_file_types = ['refpix']

    def process(self, input):

        with models.open(input) as input_model:
            if input_model.meta.exposure.readpatt is not None and \
               input_model.meta.exposure.readpatt.find("IRS2") >= 0:
                self.irs2_name = self.get_reference_file(input_model, 'irs2')
                self.log.info('Using IRS2 reference file: %s' % self.irs2_name)

                # Check for a valid reference file
                if self.irs2_name == 'N/A':
                    self.log.warning('No IRS2 reference file found')
                    self.log.warning('RefPix step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.refpix = 'SKIPPED'
                    input_model.close()
                    return result

                irs2_model = models.IRS2Model(self.irs2_name)
                result = irs2_subtract_reference.correct_model(input_model,
                                                               irs2_model)
                irs2_model.close()
            else:
                self.log.info('use_side_ref_pixels = %s' %
                              (self.use_side_ref_pixels,))
                self.log.info('odd_even_columns = %s' %
                              (self.odd_even_columns,))
                self.log.info('side_smoothing_length = %d' %
                              (self.side_smoothing_length,))
                self.log.info('side_gain = %f' % (self.side_gain,))
                self.log.info('odd_even_rows = %s' % (self.odd_even_rows,))
                result = reference_pixels.correct_model(input_model,
                                                    self.odd_even_columns,
                                                    self.use_side_ref_pixels,
                                                    self.side_smoothing_length,
                                                    self.side_gain,
                                                    self.odd_even_rows)

        result.meta.cal_step.refpix = 'COMPLETE' 

        return result

def refpix_correction(input):
    a = RefPixStep()
    result = a.process(input)
    return result
