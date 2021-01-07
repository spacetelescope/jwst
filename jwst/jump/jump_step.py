#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from .jump import detect_jumps
import time

__all__ = ["JumpStep"]


class JumpStep(Step):
    """
    JumpStep: Performs CR/jump detection on each ramp integration within an
    exposure. The 2-point difference method is applied.
    """

    spec = """
        rejection_threshold = float(default=4.0,min=0) # CR sigma rejection threshold
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
        flag_4_neighbors = boolean(default=True) # flag the four perpendicular neighbors of each CR
        max_jump_to_flag_neighbors = float(default=200) # maximum jump sigma that will trigger neighbor flagging
        min_jump_to_flag_neighbors = float(default=10) # minimum jump sigma that will trigger neighbor flagging
    """

    reference_file_types = ['gain', 'readnoise']

    def process(self, input):

        with datamodels.RampModel(input) as input_model:
            tstart = time.time()
            # Check for an input model with NGROUPS<=2
            ngroups = input_model.data.shape[1]
            if ngroups <= 4:
                self.log.warning('Can not apply jump detection when NGROUPS<=4;')
                self.log.warning('Jump step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.jump = 'SKIPPED'
                return result

            # Retrieve the parameter values
            rej_thresh = self.rejection_threshold
            max_cores = self.maximum_cores
            max_jump_to_flag_neighbors = self.max_jump_to_flag_neighbors
            min_jump_to_flag_neighbors = self.min_jump_to_flag_neighbors
            flag_4_neighbors = self.flag_4_neighbors

            self.log.info('CR rejection threshold = %g sigma', rej_thresh)
            if self.maximum_cores != 'none':
                self.log.info('Maximum cores to use = %s', max_cores)

            # Get the gain and readnoise reference files
            gain_filename = self.get_reference_file(input_model, 'gain')
            self.log.info('Using GAIN reference file: %s', gain_filename)

            gain_model = datamodels.GainModel( gain_filename )

            readnoise_filename = self.get_reference_file(input_model,
                                                          'readnoise')
            self.log.info('Using READNOISE reference file: %s',
                          readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel( readnoise_filename )

            # Call the jump detection routine
            result = detect_jumps(input_model, gain_model, readnoise_model,
                                rej_thresh, max_cores,
                                max_jump_to_flag_neighbors, min_jump_to_flag_neighbors,
                                flag_4_neighbors)

            gain_model.close()
            readnoise_model.close()
            tstop = time.time()
            self.log.info('The execution time in seconds: %f', tstop - tstart)

        result.meta.cal_step.jump = 'COMPLETE'

        return result
