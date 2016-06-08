from jwst.stpipe import Step
from jwst import datamodels
from . import rscd_sub


class RSCD_Step(Step):
    """
    RSCD_Step: Performs an RSCD correction by adding a function of time,
    frame by frame, to a copy of the input science data model.
    """

    reference_file_types = ['rscd']

    def process(self, input): 

        # Open the input data model
        with models.open(input) as input_model:

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
                    result = input_model.copy()
                    result.meta.cal_step.rscd = 'SKIPPED'
                    return result

                # Open the rscd ref file data model
                rscd_model = models.RSCD_Model(self.rscd_name)

                # Do the rscd correction subtraction
                result = rscd_sub.do_correction(input_model, rscd_model)

                # Close the reference file and update the step status
                rscd_model.close()
                result.meta.cal_step.rscd = 'COMPLETE' 

            else:
                self.log.warning('RSCD Correction is only for MIRI data')
                self.log.warning('RSCD step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.rscd = 'SKIPPED' 

        return result
