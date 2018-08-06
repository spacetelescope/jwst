from ..stpipe import Step
from .. import datamodels
from . import ami_analyze


class AmiAnalyzeStep(Step):
    """
    AmiAnalyzeStep: Performs analysis of an AMI mode exposure by
    applying the LG algorithm.
    """

    spec = """
        oversample = integer(default=3, min=1)  # Oversampling factor
        rotation = float(default=0.0)           # Rotation initial guess [deg]
    """

    reference_file_types = ['throughput']

    def process(self, input):
        """
        Performs analysis of an AMI mode exposure by applying the LG algorithm.

        Parameters
        ----------
        input: string
            input file name

        Returns
        -------
        result: AmiLgModel object
            AMI image to which the LG fringe detection has been applied
        """

        # Retrieve the parameter values
        oversample = self.oversample
        rotate = self.rotation
        self.log.info('Oversampling factor = %d', oversample)
        self.log.info('Initial rotation guess = %g deg', rotate)

        # Open the input data model
        with datamodels.ImageModel(input) as input_model:

            # Get the name of the filter throughput reference file to use
            self.filter_ref_name = self.get_reference_file(input_model,
                                   'throughput')
            self.log.info('Using filter throughput reference file %s',
                           self.filter_ref_name)

            # Check for a valid reference file
            if self.filter_ref_name == 'N/A':
                self.log.warning('No THROUGHPUT reference file found')
                self.log.warning('AMI analyze step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.ami_analyze = 'SKIPPED'
                return result

            # Open the filter throughput reference file
            filter_model = datamodels.ThroughputModel(self.filter_ref_name)

            # Do the LG analysis on the input image
            result = ami_analyze.apply_LG(input_model, filter_model,
                                          oversample, rotate)

            # Close the reference file and update the step status
            filter_model.close()
            result.meta.cal_step.ami_analyze = 'COMPLETE'

        return result
