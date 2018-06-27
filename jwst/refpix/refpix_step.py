#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import reference_pixels
from . import irs2_subtract_reference


__all__ = ["RefPixStep", "refpix_correction"]


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

        with datamodels.RampModel(input) as input_model:
            if input_model.meta.exposure.readpatt is not None and \
               input_model.meta.exposure.readpatt.find("IRS2") >= 0:
                self.irs2_name = self.get_reference_file(input_model, 'refpix')
                self.log.info('Using refpix reference file: %s' %
                              self.irs2_name)

                # Check for a valid reference file
                if self.irs2_name == 'N/A':
                    self.log.warning('No refpix reference file found')
                    self.log.warning('RefPix step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.refpix = 'SKIPPED'
                    input_model.close()
                    return result

                irs2_model = datamodels.IRS2Model(self.irs2_name)
                result = irs2_subtract_reference.correct_model(input_model,
                                                               irs2_model)
                if result.meta.cal_step.refpix != 'SKIPPED':
                    result.meta.cal_step.refpix = 'COMPLETE'
                irs2_model.close()
                return result
            else:
                self.log.info('use_side_ref_pixels = %s' %
                              (self.use_side_ref_pixels,))
                self.log.info('odd_even_columns = %s' %
                              (self.odd_even_columns,))
                self.log.info('side_smoothing_length = %d' %
                              (self.side_smoothing_length,))
                self.log.info('side_gain = %f' % (self.side_gain,))
                self.log.info('odd_even_rows = %s' % (self.odd_even_rows,))
                result = reference_pixels.correct_model(input_model.copy(),
                                                        self.odd_even_columns,
                                                        self.use_side_ref_pixels,
                                                        self.side_smoothing_length,
                                                        self.side_gain,
                                                        self.odd_even_rows)
                #
                # This returns None if there are NaNs in the output data, as would
                # be the case if the code is unable to calculate the reference value
                # due to a lack of valid reference pixels
                if result is not None:
                    if input_model.meta.subarray.name == 'FULL':
                        result.meta.cal_step.refpix = 'COMPLETE'
                    else:
                        result.meta.cal_step.refpix = 'SKIPPED'

                    return result
                else:
                    self.log.warning("Invalid reference pixels, refpix step skipped")
                    input_model.meta.cal_step.refpix = 'SKIPPED'
                    return input_model

def refpix_correction(input):
    a = RefPixStep()
    result = a.process(input)
    return result
