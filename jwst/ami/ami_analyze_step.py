from ..stpipe import Step
from .. import datamodels
from . import ami_analyze

__all__ = ["AmiAnalyzeStep"]


class AmiAnalyzeStep(Step):
    """Performs analysis of an AMI mode exposure by applying the LG algorithm.
    """

    class_alias = "ami_analyze"

    spec = """
        oversample = integer(default=3, min=1)  # Oversampling factor
        rotation = float(default=0.0)           # Rotation initial guess [deg]
        psf_offset = string(default='0.0 0.0') # Psf offset values to use to create the model array
        rotation_search = string(default='-3 3 1') # Rotation search parameters: start, stop, step
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

        # pull out parameters that are strings and change to floats
        psf_offset = [float(a) for a in self.psf_offset.split()]
        rotsearch_parameters = [float(a) for a in self.rotation_search.split()]

        self.log.info(f'Oversampling factor = {oversample}')
        self.log.info(f'Initial rotation guess = {rotate} deg')
        self.log.info(f'Initial values to use for psf offset = {psf_offset}')

        # Open the input data model
        try:
            input_model = datamodels.ImageModel(input)
        except ValueError as err:
            raise RuntimeError(f"{err}. Input must be a 2D ImageModel.")

        # check for 2D data array
        if len(input_model.data.shape) != 2:
            raise RuntimeError("Only 2D ImageModel data can be processed")

        # Get the name of the filter throughput reference file to use
        throughput_reffile = self.get_reference_file(input_model, 'throughput')
        self.log.info(f'Using filter throughput reference file {throughput_reffile}')

        # Check for a valid reference file
        if throughput_reffile == 'N/A':
            self.log.warning('No THROUGHPUT reference file found')
            self.log.warning('AMI analyze step will be skipped')
            raise RuntimeError("No throughput reference file found. "
                               "ami_analyze cannot continue.")

        # Open the filter throughput reference file
        throughput_model = datamodels.ThroughputModel(throughput_reffile)

        # Apply the LG+ methods to the data
        result = ami_analyze.apply_LG_plus(input_model, throughput_model,
                                           oversample, rotate,
                                           psf_offset,
                                           rotsearch_parameters)

        # Close the reference file and update the step status
        throughput_model.close()
        result.meta.cal_step.ami_analyze = 'COMPLETE'

        return result
