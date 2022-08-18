#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from .jump import run_detect_jumps
import time

__all__ = ["JumpStep"]


class JumpStep(Step):
    """
    JumpStep: Performs CR/jump detection on each ramp integration within an
    exposure. The 2-point difference method is applied.
    """

    spec = """
        rejection_threshold = float(default=4.0,min=0) # CR sigma rejection threshold
        three_group_rejection_threshold = float(default=6.0,min=0) # CR sigma rejection threshold
        four_group_rejection_threshold = float(default=5.0,min=0) # CR sigma rejection threshold
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
        flag_4_neighbors = boolean(default=True) # flag the four perpendicular neighbors of each CR
        max_jump_to_flag_neighbors = float(default=1000) # maximum jump sigma that will trigger neighbor flagging
        min_jump_to_flag_neighbors = float(default=10) # minimum jump sigma that will trigger neighbor flagging
        after_jump_flag_dn1 = float(default=0) # 1st flag groups after jump above DN threshold
        after_jump_flag_time1 = float(default=0) # 1st flag groups after jump groups within specified time
        after_jump_flag_dn2 = float(default=0) # 2nd flag groups after jump above DN threshold
        after_jump_flag_time2 = float(default=0) # 2nd flag groups after jump groups within specified time
    """

    reference_file_types = ['gain', 'readnoise']

    class_alias = 'jump'

    def process(self, input):

        with datamodels.RampModel(input) as input_model:
            tstart = time.time()
            # Check for an input model with NGROUPS<=2
            ngroups = input_model.data.shape[1]
            if ngroups <= 2:
                self.log.warning('Cannot apply jump detection when NGROUPS<=2;')
                self.log.warning('Jump step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.jump = 'SKIPPED'
                return result

            # Retrieve the parameter values
            rej_thresh = self.rejection_threshold
            three_grp_rej_thresh = self.three_group_rejection_threshold
            four_grp_rej_thresh = self.four_group_rejection_threshold
            max_cores = self.maximum_cores
            max_jump_to_flag_neighbors = self.max_jump_to_flag_neighbors
            min_jump_to_flag_neighbors = self.min_jump_to_flag_neighbors
            flag_4_neighbors = self.flag_4_neighbors
            after_jump_flag_dn1 = self.after_jump_flag_dn1
            after_jump_flag_time1 = self.after_jump_flag_time1
            after_jump_flag_dn2 = self.after_jump_flag_dn2
            after_jump_flag_time2 = self.after_jump_flag_time2

            self.log.info('CR rejection threshold = %g sigma', rej_thresh)
            if self.maximum_cores != 'none':
                self.log.info('Maximum cores to use = %s', max_cores)

            # Get the gain and readnoise reference files
            gain_filename = self.get_reference_file(input_model, 'gain')
            self.log.info('Using GAIN reference file: %s', gain_filename)

            gain_model = datamodels.GainModel(gain_filename)

            readnoise_filename = self.get_reference_file(input_model,
                                                         'readnoise')
            self.log.info('Using READNOISE reference file: %s',
                          readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)

            # Call the jump detection routine
            result = run_detect_jumps(input_model, gain_model, readnoise_model,
                                      rej_thresh, three_grp_rej_thresh, four_grp_rej_thresh, max_cores,
                                      max_jump_to_flag_neighbors, min_jump_to_flag_neighbors,
                                      flag_4_neighbors,
                                      after_jump_flag_dn1,
                                      after_jump_flag_time1,
                                      after_jump_flag_dn2,
                                      after_jump_flag_time2)
            gain_model.close()
            readnoise_model.close()
            tstop = time.time()
            self.log.info('The execution time in seconds: %f', tstop - tstart)

        result.meta.cal_step.jump = 'COMPLETE'

        return result
