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
        """ # noqa: E501

    #  TBD - only do this for the 2nd+ integrations
    #  do nothing for single integration exposures

    reference_file_types = ['rscd']

    def process(self, step_input):

        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:

            # check the data is MIRI data
            detector = input_model.meta.instrument.detector
            if not detector.startswith('MIR'):
                self.log.warning('RSCD correction is only for MIRI data')
                self.log.warning('RSCD step will be skipped')
                input_model.meta.cal_step.rscd = 'SKIPPED'
                return input_model

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

            # Work on a copy
            result = input_model.copy()

            # Do the rscd correction
            result = rscd_sub.do_correction(result, rscd_model, self.type)

            # Cleanup
            del rscd_model

        return result
