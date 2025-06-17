#! /usr/bin/env python

from jwst.stpipe import Step
import time

import numpy as np

from jwst.lib import reffile_utils

from stcal.jump.jump_class import JumpData
from stcal.jump.jump import detect_jumps_data

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

__all__ = ["JumpStep"]


class JumpStep(Step):
    """Step class to perform just detection using two point difference."""

    spec = """
        rejection_threshold = float(default=4.0,min=0) # CR sigma rejection threshold
        three_group_rejection_threshold = float(default=6.0,min=0) # CR sigma rejection threshold
        four_group_rejection_threshold = float(default=5.0,min=0) # CR sigma rejection threshold
        maximum_cores = string(default='1') # cores for multiprocessing. Can be an integer, 'half', 'quarter', or 'all'
        flag_4_neighbors = boolean(default=True) # flag the four perpendicular neighbors of each CR
        max_jump_to_flag_neighbors = float(default=1000) # maximum jump sigma that will trigger neighbor flagging
        min_jump_to_flag_neighbors = float(default=10) # minimum jump sigma that will trigger neighbor flagging
        after_jump_flag_dn1 = float(default=0) # 1st flag groups after jump above DN threshold
        after_jump_flag_time1 = float(default=0) # 1st flag groups after jump groups within specified time
        after_jump_flag_dn2 = float(default=0) # 2nd flag groups after jump above DN threshold
        after_jump_flag_time2 = float(default=0) # 2nd flag groups after jump groups within specified time
        expand_large_events = boolean(default=False) # Turns on Snowball detector for NIR detectors
        min_sat_area = float(default=1.0) # minimum required area for the central saturation of snowballs
        min_jump_area = float(default=5.0) # minimum area to trigger large events processing
        expand_factor = float(default=2.0) # The expansion factor for the enclosing circles or ellipses
        use_ellipses = boolean(default=False) # deprecated
        sat_required_snowball = boolean(default=True) # Require the center of snowballs to be saturated
        min_sat_radius_extend = float(default=2.5) # The min radius of the sat core to trigger the extension of the core
        sat_expand = integer(default=2) # Number of pixels to add to the radius of the saturated core of snowballs
        edge_size = integer(default=25) # Distance from detector edge where a saturated core is not required for snowball detection
        mask_snowball_core_next_int = boolean(default=True) # Flag saturated cores of snowballs in the next integration?
        snowball_time_masked_next_int = integer(default=4000) # Time in seconds over which saturated cores are flagged in next integration
        find_showers = boolean(default=False) # Apply MIRI shower flagging?
        max_shower_amplitude = float(default=4) # Maximum MIRI shower amplitude in DN/s
        extend_snr_threshold = float(default=1.2) # The SNR minimum for detection of extended showers in MIRI
        extend_min_area = integer(default=90) # Min area of emission after convolution for the detection of showers
        extend_inner_radius = float(default=1) # Inner radius of the ring_2D_kernel used for convolution
        extend_outer_radius = float(default=2.6) # Outer radius of the ring_2D_Kernel used for convolution
        extend_ellipse_expand_ratio = float(default=1.1) # Expand the radius of the ellipse fit to the extended emission
        time_masked_after_shower = float(default=15) # Seconds to flag as jump after a detected extended emission
        min_diffs_single_pass = integer(default=10) # The minimum number of differences needed to skip the iterative flagging of jumps.
        max_extended_radius = integer(default=200) # The maximum radius of an extended snowball or shower
        minimum_groups = integer(default=3) # The minimum number of groups to perform jump detection using sigma clipping
        minimum_sigclip_groups = integer(default=100) # The minimum number of groups to switch to sigma clipping
        only_use_ints = boolean(default=True) # In sigclip only compare the same group across ints, if False compare all groups
    """  # noqa: E501

    reference_file_types = ["gain", "readnoise"]

    class_alias = "jump"

    def process(self, step_input):
        """
        Step method to execute step computations.

        Parameters
        ----------
        step_input : RampModel
            The ramp model input from the previous step.

        Returns
        -------
        result : RampModel
            The ramp model with jump step as COMPLETE and jumps detected or
            the jump step is SKIPPED.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            tstart = time.time()

            # Check for an input model with NGROUPS<=2
            nints, ngroups, nrows, ncols = input_model.data.shape
            if ngroups <= 2:
                self.log.warning("Cannot apply jump detection when NGROUPS<=2;")
                self.log.warning("Jump step will be skipped")
                input_model.meta.cal_step.jump = "SKIPPED"
                return input_model

            self.log.info("CR rejection threshold = %g sigma", self.rejection_threshold)
            if self.maximum_cores != "none":
                self.log.info("Maximum cores to use = %s", self.maximum_cores)

            # Detect jumps using a copy of the input data model.
            result = input_model.copy()
            jump_data = self._setup_jump_data(result)
            new_gdq, new_pdq, number_crs, number_extended_events, stddev = detect_jumps_data(
                jump_data
            )

            # Update the DQ arrays of the output model with the jump detection results
            result.groupdq = new_gdq
            result.pixeldq = new_pdq

            # determine the number of groups with all pixels set to DO_NOT_USE
            dnu_flag = dqflags.pixel["DO_NOT_USE"]
            num_flagged_grps = 0
            for integ in range(nints):
                for grp in range(ngroups):
                    if np.all(np.bitwise_and(result.groupdq[integ, grp, :, :], dnu_flag)):
                        num_flagged_grps += 1

            total_groups = nints * ngroups - num_flagged_grps - nints
            if total_groups >= 1:
                total_time = result.meta.exposure.group_time * total_groups
                total_pixels = nrows * ncols

                crs = 1000 * number_crs / (total_time * total_pixels)
                result.meta.exposure.primary_cosmic_rays = crs

                events = 1e6 * number_extended_events / (total_time * total_pixels)
                result.meta.exposure.extended_emission_events = events

            tstop = time.time()
            self.log.info("The execution time in seconds: %f", tstop - tstart)

            result.meta.cal_step.jump = "COMPLETE"

        return result

    def _setup_jump_data(self, result):
        """
        Create a JumpData instance to be used by STCAL jump.

        Parameters
        ----------
        result : RampModel
            The ramp model input from the previous step.

        Returns
        -------
        jump_data : JumpData
            The data container to be used to run the STCAL detect_jumps_data.
        """
        # Get the gain and readnoise reference files
        gain_filename = self.get_reference_file(result, "gain")
        self.log.info("Using GAIN reference file: %s", gain_filename)
        readnoise_filename = self.get_reference_file(result, "readnoise")
        self.log.info("Using READNOISE reference file: %s", readnoise_filename)

        with (
            datamodels.ReadnoiseModel(readnoise_filename) as rnoise_m,
            datamodels.GainModel(gain_filename) as gain_m,
        ):
            # Get 2D gain and read noise values from their respective models
            if reffile_utils.ref_matches_sci(result, gain_m):
                gain_2d = gain_m.data
            else:
                self.log.info("Extracting gain subarray to match science data")
                gain_2d = reffile_utils.get_subarray_model(result, gain_m).data

            if reffile_utils.ref_matches_sci(result, rnoise_m):
                rnoise_2d = rnoise_m.data
            else:
                self.log.info("Extracting readnoise subarray to match science data")
                rnoise_2d = reffile_utils.get_subarray_model(result, rnoise_m).data

        # Instantiate a JumpData class and populate it based on the input RampModel.
        jump_data = JumpData(result, gain_2d, rnoise_2d, dqflags.pixel)

        jump_data.set_detection_settings(
            self.rejection_threshold,
            self.three_group_rejection_threshold,
            self.four_group_rejection_threshold,
            self.max_jump_to_flag_neighbors,
            self.min_jump_to_flag_neighbors,
            self.flag_4_neighbors,
        )

        # determine the number of groups that correspond to the after_jump times
        # needed because the group time is not passed to detect_jumps_data
        gtime = result.meta.exposure.group_time
        after_jump_flag_n1 = int(self.after_jump_flag_time1 // gtime)
        after_jump_flag_n2 = int(self.after_jump_flag_time2 // gtime)

        jump_data.set_after_jump(
            self.after_jump_flag_dn1,
            after_jump_flag_n1,
            self.after_jump_flag_dn2,
            after_jump_flag_n2,
        )

        sat_expand = self.sat_expand * 2
        jump_data.set_snowball_info(
            self.expand_large_events,
            self.min_jump_area,
            self.min_sat_area,
            self.expand_factor,
            self.sat_required_snowball,
            self.min_sat_radius_extend,
            sat_expand,
            self.edge_size,
        )

        max_extended_radius = self.max_extended_radius * 2
        jump_data.set_shower_info(
            self.find_showers,
            self.extend_snr_threshold,
            self.extend_min_area,
            self.extend_inner_radius,
            self.extend_outer_radius,
            self.extend_ellipse_expand_ratio,
            self.min_diffs_single_pass,
            max_extended_radius,
        )

        jump_data.set_sigma_clipping_info(
            self.minimum_groups, self.minimum_sigclip_groups, self.only_use_ints
        )

        jump_data.max_cores = self.maximum_cores
        jump_data.grps_masked_after_shower = int(self.time_masked_after_shower // gtime)
        jump_data.mask_persist_grps_next_int = self.mask_snowball_core_next_int
        jump_data.persist_grps_flagged = int(self.snowball_time_masked_next_int // gtime)
        jump_data.max_shower_amplitude = jump_data.max_shower_amplitude * gtime

        return jump_data
