#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
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
        min_sat_area = float(default=1.0) # minimum required area for the central saturation of snowballs
        min_jump_area = float(default=5.0) # minimum area to trigger large events processing
        expand_factor = float(default=2.0) # The expansion factor for the enclosing circles or ellipses
        use_ellipses = boolean(default=False) # deprecated
        sat_required_snowball = boolean(default=True) # Require the center of snowballs to be saturated
        min_sat_radius_extend = float(default=2.5) # The min radius of the sat core to trigger the extension of the core
        sat_expand = integer(default=2) Number of pixels to add to the radius of the saturated core of snowballs
        expand_large_events = boolean(default=False) # Turns on Snowball detector for NIR detectors
        find_showers = boolean(default=False) Turn on shower flagging for MIRI
        edge_size = integer(default=25) # Size of region on the edges of NIR detectors where a sat core is not required
        extend_snr_threshold = float(default=1.2) The SNR minimum for detection of extended showers in MIRI
        extend_min_area = integer(default=90) Min area of emission after convolution for the detection of showers
        extend_inner_radius = float(default=1) Inner radius of the ring_2D_kernel used for convolution
        extend_outer_radius = float(default=2.6) Outer radius of the ring_2D_Kernel used for convolution
        extend_ellipse_expand_ratio = float(default=1.1) Expand the radius of the ellipse fit to the extended emission
        time_masked_after_shower = float(default=15) Seconds to flag as jump after a detected extended emission
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
            min_sat_area = self.min_sat_area
            min_jump_area = self.min_jump_area
            expand_factor = self.expand_factor
            use_ellipses = self.use_ellipses
            sat_required_snowball = self.sat_required_snowball
            expand_large_events = self.expand_large_events
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
                                      after_jump_flag_time2,
                                      min_sat_area=min_sat_area, min_jump_area=min_jump_area,
                                      expand_factor=expand_factor, use_ellipses=use_ellipses,
                                      min_sat_radius_extend=self.min_sat_radius_extend,
                                      sat_required_snowball=sat_required_snowball, sat_expand=self.sat_expand,
                                      expand_large_events=expand_large_events, find_showers=self.find_showers,
                                      edge_size=self.edge_size, extend_snr_threshold=self.extend_snr_threshold,
                                      extend_min_area=self.extend_min_area,
                                      extend_inner_radius=self.extend_inner_radius,
                                      extend_outer_radius=self.extend_outer_radius,
                                      extend_ellipse_expand_ratio=self.extend_ellipse_expand_ratio,
                                      time_masked_after_shower=self.time_masked_after_shower)

            gain_model.close()
            readnoise_model.close()
            tstop = time.time()
            self.log.info('The execution time in seconds: %f', tstop - tstart)

        result.meta.cal_step.jump = 'COMPLETE'

        return result
