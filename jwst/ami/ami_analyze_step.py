from ..stpipe import Step
from .. import datamodels
from . import ami_analyze

__all__ = ["AmiAnalyzeStep"]


class AmiAnalyzeStep(Step):
    """Performs analysis of an AMI mode exposure by applying the LG algorithm.
    """
    spec = """
        oversample = integer(default=3, min=1)  # Oversampling factor
        rotation = float(default=0.0)           # Rotation initial guess [deg]
        set_psf_offset_find_rotation = string(default='0.0 0.0') # Psf offset values to use to create the model array
        set_rotsearch_d = string(default='-3 3.1 1.') # Set of values to use for the rotation search
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

        psf_offset = self.set_psf_offset_find_rotation.split(' ')
        psf_offset_find_rotation = (float(psf_offset[0]), float(psf_offset[1]))
        rotsearch = self.set_rotsearch_d.split(' ')
        rotsearch_parameters= [float(rotsearch[0]), float(rotsearch[1]), float(rotsearch[2])]

        self.log.info(f'Oversampling factor =  {oversample}')
        self.log.info(f'Initial rotation guess = {rotate} deg')
        self.log.info(f'Initial values to use for psf offset = {psf_offset_find_rotation}')
        self.log.info(f'Initial values to use for rotation search {rotsearch_parameters}')

        # Open the input data model
        try:
            with datamodels.ImageModel(input) as input_model:
                if len(input_model.data.shape) != 2:
                    raise RuntimeError("Only 2D ImageModel data can be processed.")
                # Get the name of the filter throughput reference file to use
                throughput_reffile = self.get_reference_file(input_model,
                                       'throughput')

                self.log.info('Using filter throughput reference file %s',
                               throughput_reffile)

                # Check for a valid reference file
                if throughput_reffile == 'N/A':
                    self.log.warning('No THROUGHPUT reference file found')
                    self.log.warning('AMI analyze step will be skipped')
                    raise RuntimeError("No throughput reference file found. "
                                       "ami_analyze cannot continue.")

                # Open the filter throughput reference file
                throughput_model = datamodels.ThroughputModel(throughput_reffile)

                result = ami_analyze.apply_LG_plus(input_model, throughput_model,
                                                   oversample, rotate,
                                                   psf_offset_find_rotation,
                                                   rotsearch_parameters)

                # Close the reference file and update the step status
                throughput_model.close()
                result.meta.cal_step.ami_analyze = 'COMPLETE'

            return result

        # If _calints CubeModel input, handle as RuntimeError
        except ValueError as err:
            raise RuntimeError("Only 2D ImageModel data can be processed.") from err
