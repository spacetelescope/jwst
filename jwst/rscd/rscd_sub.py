#
#  Module for the RSCD correction for MIRI science data
#

import logging

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)

__all__ = ["do_correction", "correction_skip_groups", "get_rscd_parameters"]


def do_correction(output_model, rscd_model, bright_use_group1=False):
    """
    Set the initial groups of an integration for  MIRI data to 'DO_NOT_USE'.

    The number of initial groups to set to 'DO_NOT_USE' is read in from the RSCD reference
    file. The number of groups to skip is integration dependent. The first integration has
    a value defined in the reference file and the second and higher integrations have a
    separate value in the reference file.

    Parameters
    ----------
    output_model : RampModel
        Input ramp datamodel

    rscd_model : RSCDModel
        RSCD reference datamodel

    bright_use_group1 : bool
        If True, setting first group data quality flag to DO_NOT_USE will not be
        done for pixels that have the saturation flag set for group 3.

    Returns
    -------
    output_model : RampModel
        Ramp datamodel with RSCD affected groups flagged as DO_NOT_USE
    """
    # Retrieve the reference parameters for this exposure type
    param = get_rscd_parameters(output_model, rscd_model)

    if not bool(param):  # empty dictionary
        log.warning(
            "READPATT, SUBARRAY combination not found in ref file: RSCD correction will be skipped"
        )
        output_model.meta.cal_step.rscd = "SKIPPED"
        return output_model

    group_skip_int1 = param["skip_int1"]  # integration 1
    group_skip_int2p = param["skip_int2p"]  # integration 2,  plus higher integrations
    log.info(f" # groups from RSCD reference file for int 1 to flag  {group_skip_int1}")
    log.info(f" # groups from RSCD reference file for int 2 and higher to flag  {group_skip_int2p}")
    output_model = correction_skip_groups(
        output_model, group_skip_int1, group_skip_int2p, bright_use_group1
    )
    return output_model


def correction_skip_groups(output, group_skip_int1, group_skip_int2p, bright_use_group1=False):
    """
    Set the initial groups in integration to DO_NOT_USE to skip groups affected by RSCD effect.

    Parameters
    ----------
    output : RampModel
        Science data to be flagged

    group_skip_int1 : int
        Number of groups to skip at the beginning of the ramp for integration 1

    group_skip_int2p : int
        Number of groups to skip at the beginning of the ramp for integration 2 and higher

    Returns
    -------
    output_model: RampModel
        Ramp datamodel with RSCD affected groups flagged as DO_NOT_USE
    """
    # General exposure parameters
    sci_ngroups = output.meta.exposure.ngroups
    sci_nints = output.meta.exposure.nints

    # values defined for segmented data
    sci_int_start = output.meta.exposure.integration_start

    if sci_int_start is None:
        sci_int_start = 1

    log.debug(f"RSCD correction using: nints={sci_nints}, ngroups={sci_ngroups}")
    log.debug(f"The first integration in the data is integration: {sci_int_start}")

    # We need at least 3 groups for the Bright Use Group1 (non-flagging)
    # After an RSCD correction we have to have 3 usable groups.
    # To flag groups affected by the RSCD correction we need at least 4 groups and that would allow
    # us to only flag 1 group. If there are only 3 groups (which should not occur for MIRI data)
    # we would need to skip flagging RSCD affected groups.
    if sci_ngroups < 3:
        log.warning("Too few groups to apply RSCD correction")
        log.warning("RSCD step will be skipped")
        output.meta.cal_step.rscd = "SKIPPED"
        return output

    # Will we have at least 3 groups left after skipping RSCD affected groups ? If not reduced
    # the number of groups skipped to ensure we have at least 3 groups. We need 3 because the
    # last frame correction will flag the last group as DO_NOT_USE.

    # checks for integration 1:
    if sci_ngroups < (group_skip_int1 + 3):
        max_groups_skip = max(0, sci_ngroups - 3)

        if max_groups_skip != group_skip_int1:
            log.info(f"Changing the # of groups to skip in int 1 to {max_groups_skip}")
            group_skip_int1 = max_groups_skip

    # checks for integration 2
    if sci_nints > 1:
        if sci_ngroups < (group_skip_int2p + 3):  # second and higher checks on the number of groups
            max_groups_skip = max(0, sci_ngroups - 3)
            if max_groups_skip != group_skip_int2p:
                group_skip_int2p = max_groups_skip
                log.info(
                    f"Changing the # of groups to skip in int 2 and higher to {max_groups_skip}"
                )

    # For segmented data the first integration in the file may not be the first
    # integration in the exposure. The value in meta.exposure.integration_start
    # holds the value of the first integration in the file.
    # If a correction is to be done and if ngroups > group_skip+3, then set all
    # of the GROUPDQ flags from groups 0 to group_skip to 'DO_NOT_USE'.

    # drop groups in integration 1
    # first integration
    log.info(f"Number of groups to skip for integration 1: {group_skip_int1}")

    if sci_int_start == 1:  # we have segmented data and the data starts with the first integration.
        if bright_use_group1:
            # --- Logic for retaining Group 1 in bright, saturated pixels ---

            # 1. Identify pixels where saturation happens in Group 3 (index 2).
            # In this specific case, the signal in group2-group1 is likely dominant
            # over the first frame effect, so we keep group 1.
            is_saturated_in_group3 = (output.groupdq[0, 2, :, :] & dqflags.group["SATURATED"]) > 0
            is_not_saturated_in_group1 = (
                output.groupdq[0, 0, :, :] & dqflags.group["SATURATED"]
            ) == 0

            is_bright_use_group1 = is_saturated_in_group3 & is_not_saturated_in_group1
            row_indices, col_indices = np.where(~is_bright_use_group1)
            # Following in just for testing - will remove  from code for final PR

            # Find the indices where the array is NOT saturated (where the inverted array is True)

            # only flag pixels that do not meet requirement that group3 saturats but group2 does not
            output.groupdq[0, 0:group_skip_int1, row_indices, col_indices] = np.bitwise_or(
                output.groupdq[0, 0:group_skip_int1, row_indices, col_indices],
                dqflags.group["DO_NOT_USE"],
            )

            # also set a FLUX_ESTIMATED flag in the 2D pixeldq image
            # for any group1 pixels that were kept
            output.pixeldq[is_bright_use_group1] |= dqflags.pixel["FLUX_ESTIMATED"]

            any_good_kept = (
                output.pixeldq[is_bright_use_group1] & dqflags.pixel["DO_NOT_USE"]
            ) == 0

            log.info(
                "Number of usable bright pixels with first group "
                f"not set to DO_NOT_USE: {np.sum(any_good_kept)}"
            )
            output.meta.rscd.keep_bright_firstgroup_int1 = np.sum(any_good_kept)
        else:  # RSCD flagging of every pixel
            output.groupdq[0, 0:group_skip_int1, :, :] = np.bitwise_or(
                output.groupdq[0, 0:group_skip_int1, :, :], dqflags.group["DO_NOT_USE"]
            )

    output.meta.rscd.ngroups_skip_int1 = group_skip_int1

    # Now are we flagging the second and higher integration ?
    int_start = 1
    if sci_int_start != 1:  # we have segmented data and we are not on the first integration
        int_start = 0

    if sci_nints > 1:
        log.info(f"Number of groups to skip for integrations 2 and higher: {group_skip_int2p}")

        if bright_use_group1:
            # --- Logic for retaining Group 1 in bright, saturated pixels ---

            # 1. Identify pixels where saturation happens in Group 3 (index 2).
            # In this specific case, the signal in group2-group1 is likely dominant
            # over the first frame effect, so we keep group 1.
            is_saturated_in_group3 = (
                output.groupdq[int_start:, 2, :, :] & dqflags.group["SATURATED"]
            ) > 0
            is_not_saturated_in_group1 = (
                output.groupdq[int_start:, 0, :, :] & dqflags.group["SATURATED"]
            ) == 0

            is_bright_use_group1 = is_saturated_in_group3 & is_not_saturated_in_group1
            # Find the indices where the array is NOT saturated (where the inverted array is True)
            int_indices, row_indices, col_indices = np.where(~is_bright_use_group1)
            int_indices = int_indices + int_start

            output.groupdq[int_indices, 0:group_skip_int2p, row_indices, col_indices] = (
                np.bitwise_or(
                    output.groupdq[int_indices, 0:group_skip_int2p, row_indices, col_indices],
                    dqflags.group["DO_NOT_USE"],
                )
            )
            # also set a FLUX_ESTIMATED flag in the 2D pixeldq image
            # for any group1 pixels that were kept
            # if we have 3 D data (more than 1 integration of data contained in is_bright_use_group1
            # need to convert to 2D data

            any_kept = np.any(is_bright_use_group1, axis=0)
            output.pixeldq[any_kept] |= dqflags.pixel["FLUX_ESTIMATED"]

            any_good_kept = (output.pixeldq[any_kept] & dqflags.pixel["DO_NOT_USE"]) == 0
            log.info(
                "Number of usable bright pixels with first group "
                f"not set to DO_NOT_USE: {np.sum(any_good_kept)}"
            )
            output.meta.rscd.keep_bright_firstgroup_int2p = np.sum(any_good_kept)

        else:
            output.groupdq[int_start:, 0:group_skip_int2p, :, :] = np.bitwise_or(
                output.groupdq[int_start:, 0:group_skip_int2p, :, :], dqflags.group["DO_NOT_USE"]
            )

        output.meta.rscd.ngroups_skip_int2p = group_skip_int2p
        output.meta.cal_step.rscd = "COMPLETE"

    return output


def get_rscd_parameters(input_model, rscd_model):
    """
    Read in the parameters from the reference file and store the parameters in a dictionary.

    Parameters
    ----------
    input_model : RampModel
        Science data to be flagged

    rscd_model : RSCDModel
        RSCD reference file data

    Returns
    -------
    param : dict
        Dictionary of parameters
    """
    # Reference file parameters held in dictionary: param
    param = {}

    # read in the type of data from the input model (FAST,SLOW,FULL,SUBARRAY)
    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name

    # Check for old values of the MIRI LRS slitless subarray name
    # in the science data and change to the new
    if subarray.upper() == "SUBPRISM":
        subarray = "SLITLESSPRISM"

    # read table 1: containing the number of groups to skip
    for tabdata in rscd_model.rscd_group_skip_table:
        subarray_table = tabdata["subarray"]
        readpatt_table = tabdata["readpatt"]
        group_skip_table_int2p = tabdata["group_skip"]  # integration 2 and higher (+)
        group_skip_table_int1 = tabdata["group_skip1"]

        if subarray_table == subarray and readpatt_table == readpatt:
            param["skip_int1"] = group_skip_table_int1
            param["skip_int2p"] = group_skip_table_int2p  # integration 2 and higher
            break

    return param
