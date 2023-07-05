import logging
import numpy as np
from stcal.jump.jump import detect_jumps
from stdatamodels.jwst.datamodels import dqflags

from ..lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def run_detect_jumps(input_model, gain_model, readnoise_model,
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
                     time_masked_after_shower=30,
                     max_extended_radius=200,
                     minimum_groups=3,
                     minimum_sigclip_groups=100,
                     only_use_ints=True
                     ):

    # Runs `detect_jumps` in stcal

    # extract data and info from input_model to pass to detect_jumps
    frames_per_group = input_model.meta.exposure.nframes
    data = input_model.data
    gdq = input_model.groupdq
    pdq = input_model.pixeldq
    err = input_model.err
    output_model = input_model.copy()

    # determine the number of groups that correspond to the after_jump times
    # needed because the group time is not passed to detect_jumps
    gtime = input_model.meta.exposure.group_time
    after_jump_flag_n1 = int(after_jump_flag_time1 // gtime)
    after_jump_flag_n2 = int(after_jump_flag_time2 // gtime)
    grps_masked_after_shower = int(time_masked_after_shower // gtime)

    # Get 2D gain and read noise values from their respective models
    if reffile_utils.ref_matches_sci(input_model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info('Extracting gain subarray to match science data')
        gain_2d = reffile_utils.get_subarray_data(input_model, gain_model)

    if reffile_utils.ref_matches_sci(input_model, readnoise_model):
        readnoise_2d = readnoise_model.data
    else:
        log.info('Extracting readnoise subarray to match science data')
        readnoise_2d = reffile_utils.get_subarray_data(input_model,
                                                       readnoise_model)
    new_gdq, new_pdq, number_crs, number_extended_events, stddev\
        = detect_jumps(frames_per_group, data, gdq, pdq, err,
                                    gain_2d, readnoise_2d,
                                    rejection_thresh, three_grp_thresh,
                                    four_grp_thresh, max_cores,
                                    max_jump_to_flag_neighbors,
                                    min_jump_to_flag_neighbors,
                                    flag_4_neighbors, dqflags.pixel,
                                    after_jump_flag_dn1,
                                    after_jump_flag_n1,
                                    after_jump_flag_dn2,
                                    after_jump_flag_n2,
                                    min_sat_area=min_sat_area, min_jump_area=min_jump_area,
                                    expand_factor=expand_factor, use_ellipses=use_ellipses,
                                    min_sat_radius_extend=min_sat_radius_extend,
                                    sat_required_snowball=sat_required_snowball, sat_expand=sat_expand,
                                    expand_large_events=expand_large_events, find_showers=find_showers,
                                    edge_size=edge_size, extend_snr_threshold=extend_snr_threshold,
                                    extend_min_area=extend_min_area, extend_inner_radius=extend_inner_radius,
                                    extend_outer_radius=extend_outer_radius,
                                    extend_ellipse_expand_ratio=extend_ellipse_expand_ratio,
                                    grps_masked_after_shower=grps_masked_after_shower,
                                    max_extended_radius=max_extended_radius,
                                    minimum_groups=minimum_groups,
                                    minimum_sigclip_groups=minimum_sigclip_groups,
                                    only_use_ints=only_use_ints
                                    )


    # Update the DQ arrays of the output model with the jump detection results
    output_model.groupdq = new_gdq
    output_model.pixeldq = new_pdq
    # determine the number of groups with all pixels set to DO_NOT_USE
    dnu_flag = 1
    num_flagged_grps = 0
    for integ in range(data.shape[0]):
        for grp in range(data.shape[1]):
            if np.all(np.bitwise_and(gdq[integ, grp, :, :], dnu_flag)):
                num_flagged_grps += 1
    total_groups = data.shape[0] * data.shape[1] - num_flagged_grps - data.shape[0]
    total_time = output_model.meta.exposure.group_time * total_groups
    total_pixels = data.shape[2] * data.shape[3]
    output_model.meta.exposure.primary_cosmic_rays = 1000 * number_crs / (total_time * total_pixels)
    output_model.meta.exposure.extended_emission_events = 1e6 * number_extended_events /\
                                                         (total_time * total_pixels)

    return output_model
