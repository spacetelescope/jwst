#! /usr/bin/env python

from ..stpipe import Step
from ..lib import pipe_utils
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
            # shape[-2] will be 3200 for IRS2 data.
            if pipe_utils.is_irs2(input_model):
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
                datamodel = input_model.copy()
                status = reference_pixels.correct_model(datamodel,
                                                        self.odd_even_columns,
                                                        self.use_side_ref_pixels,
                                                        self.side_smoothing_length,
                                                        self.side_gain,
                                                        self.odd_even_rows)
                if status == reference_pixels.REFPIX_OK:
                    datamodel.meta.cal_step.refpix = 'COMPLETE'
                elif status == reference_pixels.SUBARRAY_DOESNTFIT:
                    self.log.warning("Subarray doesn't fit in full-sized array")
                    datamodel.meta.cal_step.refpix = 'SKIPPED'
                elif status == reference_pixels.BAD_REFERENCE_PIXELS:
                    self.log.warning("No valid reference pixels, refpix step skipped")
                    datamodel.meta.cal_step.refpix = 'SKIPPED'
                elif status == reference_pixels.SUBARRAY_SKIPPED:
                    datamodel.meta.cal_step.refpix = 'SKIPPED'
                return datamodel

def refpix_correction(input):
    a = RefPixStep()
    result = a.process(input)
    return result
