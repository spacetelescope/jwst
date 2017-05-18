from jwst_lib.stpipe import Step
from jwst_lib import models
from . import msaflag_open


class MSAFlagOpenStep(Step):
    """
    MSAFlagOpenStep: Flags pixels affected by MSA failed open shutters
    """

    spec = """

    """

    reference_file_types = ['badshutter']

    def process(self, input):

        # Open the input data model
        with models.open(input) as input_model:

            # Get the name of the reference file to use
            self.reference_name = self.get_reference_file(input_model,
                                                          'badshutter')
            self.log.info('Using reference file %s', self.reference_name)

            # Check for a valid reference file
            if self.reference_name == 'N/A':
                self.log.warning('No reference file found')
                self.log.warning('Step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.msaflagopen = 'SKIPPED'
                return result

            # Open the bad shutter ref file data model
            f1 = open(self.reference_name, 'r')
            shutters = json.load(f1)
            f1.close()

            # Do the DQ flagging
            result = msaflag_open.do_correction(input_model,
                                                shutters['badshutters'])

            # set the step status to complete
            result.meta.cal_step.msaflagopen = 'COMPLETE'

        return result

