from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import rscd_sub

__all__ = ["RscdStep"]


class RscdStep(Step):
    """
    RscdStep: Performs an RSCD correction to MIRI data.
    Baseline version flags the first N groups as 'DO_NOT_USE' in
    the 2nd and later integrations in a copy of the input
    science data model.
    Enhanced version is not ready nor enabled.
    """

    class_alias = "rscd"

    # allow switching between baseline and enhanced algorithms
    spec = """
         type = option('baseline','enhanced',default = 'baseline') # Type of correction
       """

    #  TBD - only do this for the 2nd+ integrations
    #  do nothing for single integration exposures

    reference_file_types = ['rscd']

    def process(self, input):

        # Open the input data model
        with datamodels.RampModel(input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if detector.startswith('MIR'):

                # Get the name of the rscd reference file to use
                self.rscd_name = self.get_reference_file(input_model, 'rscd')
                self.log.info('Using RSCD reference file %s', self.rscd_name)

                # Check for a valid reference file
                if self.rscd_name == 'N/A':
                    self.log.warning('No RSCD reference file found')
                    self.log.warning('RSCD step will be skipped')
                    input_model.meta.cal_step.rscd = 'SKIPPED'
                    return input_model

                # Load the rscd ref file data model
                rscd_model = datamodels.RSCDModel(self.rscd_name)

                # Do the rscd correction
                result = rscd_sub.do_correction(input_model, rscd_model, self.type)

                # Close the reference file
                rscd_model.close()

            else:
                self.log.warning('RSCD correction is only for MIRI data')
                self.log.warning('RSCD step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.rscd = 'SKIPPED'

        return result
