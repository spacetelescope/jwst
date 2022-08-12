import logging

from ..datamodels import dqflags
from ..lib import reffile_utils
from stcal.jump.jump import detect_jumps

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def run_detect_jumps(input_model, gain_model, readnoise_model,
                     rejection_thresh, three_grp_thresh, four_grp_thresh,
                     max_cores, max_jump_to_flag_neighbors,
                     min_jump_to_flag_neighbors, flag_4_neighbors,
                     after_jump_flag_dn1=0.0,
                     after_jump_flag_time1=0.0,
                     after_jump_flag_dn2=0.0,
                     after_jump_flag_time2=0.0):

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

    new_gdq, new_pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                                    gain_2d, readnoise_2d,
                                    rejection_thresh, three_grp_thresh,
                                    four_grp_thresh, max_cores,
                                    max_jump_to_flag_neighbors,
                                    min_jump_to_flag_neighbors,
                                    flag_4_neighbors, dqflags.pixel,
                                    after_jump_flag_dn1,
                                    after_jump_flag_n1,
                                    after_jump_flag_dn2,
                                    after_jump_flag_n2)

    # Update the DQ arrays of the output model with the jump detection results
    output_model.groupdq = new_gdq
    output_model.pixeldq = new_pdq

    return output_model
