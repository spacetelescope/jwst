# XXX
"""
The function run_detect_jumps is no longer used by the jump step.

This file is no longer needed, but it cannot simply be deleted because
there are CI tests thet call run_detect_jumps directly.  Those tests are
in test_detect_jumps.py.  The tests in that test suite should be refactored
and moved to test_jump_step.py now that run_detect_jumps is no longer used.
"""

import logging
import numpy as np
from stcal.jump.jump import detect_jumps_data
from stcal.jump.jump_class import JumpData
from stdatamodels.jwst.datamodels import dqflags

from ..lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def run_detect_jumps(
        output_model, gain_model, readnoise_model,

        rejection_thresh, three_grp_thresh, four_grp_thresh,
        max_cores, max_jump_to_flag_neighbors,
        min_jump_to_flag_neighbors, flag_4_neighbors,

        after_jump_flag_dn1=0.0,
        after_jump_flag_time1=0.0,
        after_jump_flag_dn2=0.0,
        after_jump_flag_time2=0.0,

        min_sat_area=1.0, min_jump_area=5.0, min_sat_radius_extend=2.5,
        expand_factor=2.0, use_ellipses=False,
        sat_required_snowball=True, sat_expand=2,

        expand_large_events=False, find_showers=False, edge_size=25, extend_snr_threshold=1.1,
        extend_min_area=90, extend_inner_radius=1, extend_outer_radius=2.6, extend_ellipse_expand_ratio=1.1,

        time_masked_after_shower=30, min_diffs_single_pass=10,
        max_extended_radius=200,
        minimum_groups=3,
        minimum_sigclip_groups=100,

        only_use_ints=True,
        mask_snowball_persist_next_int=True,
        snowball_time_masked_next_int=250,
        max_shower_amplitude=4
):
    # determine the number of groups that correspond to the after_jump times
    # needed because the group time is not passed to detect_jumps_data
    gtime = output_model.meta.exposure.group_time
    after_jump_flag_n1 = int(after_jump_flag_time1 // gtime)
    after_jump_flag_n2 = int(after_jump_flag_time2 // gtime)
    grps_masked_after_shower = int(time_masked_after_shower // gtime)
    snowball_grps_masked_next_int = int(snowball_time_masked_next_int // gtime)

    # Likewise, convert a max MIRI shower amplitude in DN/s to DN/group
    max_shower_amplitude = max_shower_amplitude * gtime

    # Get 2D gain and read noise values from their respective models
    if reffile_utils.ref_matches_sci(output_model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info('Extracting gain subarray to match science data')
        gain_2d = reffile_utils.get_subarray_data(output_model, gain_model)

    if reffile_utils.ref_matches_sci(output_model, readnoise_model):
        readnoise_2d = readnoise_model.data
    else:
        log.info('Extracting readnoise subarray to match science data')
        readnoise_2d = reffile_utils.get_subarray_data(output_model, readnoise_model)

    # Create and populate the RampData class
    jump_data = JumpData(output_model, gain_2d, readnoise_2d, dqflags.pixel)
    jump_data.set_detection_settings(
        rejection_thresh, three_grp_thresh, four_grp_thresh, max_jump_to_flag_neighbors,
        min_jump_to_flag_neighbors, flag_4_neighbors)
    jump_data.set_after_jump(
        after_jump_flag_dn1, after_jump_flag_n1, after_jump_flag_dn2, after_jump_flag_n2)
    jump_data.set_snowball_info(
        expand_large_events, min_jump_area, min_sat_area, expand_factor, sat_required_snowball,
        min_sat_radius_extend, sat_expand, edge_size)
    jump_data.set_shower_info(
        find_showers, extend_snr_threshold, extend_min_area, extend_inner_radius,
        extend_outer_radius, extend_ellipse_expand_ratio, min_diffs_single_pass, max_extended_radius)
    jump_data.set_sigma_clipping_info(
        minimum_groups, minimum_sigclip_groups, only_use_ints)
    jump_data.max_cores = max_cores
    jump_data.grps_masked_after_shower = int(time_masked_after_shower // gtime)
    jump_data.mask_persist_grps_next_int = mask_snowball_persist_next_int
    jump_data.persist_grps_flagged = int(snowball_time_masked_next_int // gtime)

    # Detect jumps
    new_gdq, new_pdq, number_crs, number_extended_events, stddev = detect_jumps_data(jump_data)
        
    # Update the DQ arrays of the output model with the jump detection results
    output_model.groupdq = new_gdq
    output_model.pixeldq = new_pdq

    # determine the number of groups with all pixels set to DO_NOT_USE
    dnu_flag = 1
    num_flagged_grps = 0
    datashape = np.shape(output_model.data)
    for integ in range(datashape[0]):
        for grp in range(datashape[1]):
            if np.all(np.bitwise_and(output_model.groupdq[integ, grp, :, :], dnu_flag)):
                num_flagged_grps += 1
    total_groups = datashape[0] * datashape[1] - num_flagged_grps - datashape[0]
    if total_groups >= 1:
        total_time = output_model.meta.exposure.group_time * total_groups
        total_pixels = datashape[2] * datashape[3]
        output_model.meta.exposure.primary_cosmic_rays = 1000 * number_crs / (total_time * total_pixels)
        output_model.meta.exposure.extended_emission_events = 1e6 * number_extended_events /\
                                                              (total_time * total_pixels)

    return output_model
