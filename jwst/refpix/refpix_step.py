from stdatamodels.jwst import datamodels

from ..stpipe import Step
from ..lib import pipe_utils
from . import reference_pixels
from . import irs2_subtract_reference


__all__ = ["RefPixStep"]


class RefPixStep(Step):
    """
    RefPixStep: Use reference pixels to correct bias drifts
    """

    class_alias = "refpix"

    spec = """
        odd_even_columns = boolean(default=True) # Compute reference signal separately for even/odd columns
        use_side_ref_pixels = boolean(default=True) # Use side reference pixels for reference signal for each row
        side_smoothing_length = integer(default=11) # Median window smoothing height for side reference signal
        side_gain = float(default=1.0) # Multiplicative factor for side reference signal before subtracting from rows
        odd_even_rows = boolean(default=True) # Compute reference signal separately for even- and odd-numbered rows
        ovr_corr_mitigation_ftr = float(default=3.0) # Factor to avoid overcorrection of bad reference pixels for IRS2
    """

    reference_file_types = ['refpix']

    def process(self, input):

        # Load the input science data
        with datamodels.RampModel(input) as input_model:

            if pipe_utils.is_irs2(input_model):

                # If the science data uses NIRSpec IRS2 readout mode,
                # get the necessary refpix reference file
                self.irs2_name = self.get_reference_file(input_model, 'refpix')
                self.log.info(f'Using refpix reference file: {self.irs2_name}')

                # Check for a valid reference file
                if self.irs2_name == 'N/A':
                    self.log.warning('No refpix reference file found')
                    self.log.warning('RefPix step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.refpix = 'SKIPPED'
                    input_model.close()
                    return result

                # Load the reference file into a datamodel
                irs2_model = datamodels.IRS2Model(self.irs2_name)

                # Apply the IRS2 correction scheme
                result = irs2_subtract_reference.correct_model(input_model, irs2_model,
                                  ovr_corr_mitigation_ftr=self.ovr_corr_mitigation_ftr)

                if result.meta.cal_step.refpix != 'SKIPPED':
                    result.meta.cal_step.refpix = 'COMPLETE'
                irs2_model.close()
                return result

            else:
                # Not an NRS IRS2 exposure. Do the normal refpix correction.
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
