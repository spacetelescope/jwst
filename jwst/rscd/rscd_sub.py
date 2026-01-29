#
#  Module for the RSCD correction for MIRI science data
#

import logging

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)

__all__ = [
    "do_correction",
    "correction_skip_groups",
    "get_rscd_parameters",
    "flag_rscd",
    "apply_rscd_flags",
]


def do_correction(output_model, rscd_model, bright_use_2groups=True):
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

    bright_use_2groups : bool
        If True, for saturated data adjust number of rscd groups to flag to ensure we
        have at least 2 valid groups.

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
        output_model, group_skip_int1, group_skip_int2p, bright_use_2groups
    )
    return output_model


def correction_skip_groups(output, group_skip_int1, group_skip_int2p, bright_use_2groups=True):
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

    bright_use_2groups : bool
        If True, for saturated data adjust number of rscd groups to flag to ensure we
        have at least 2 valid groups.

    Returns
    -------
    output: RampModel
        Ramp datamodel with RSCD affected groups flagged as DO_NOT_USE
    """
    # General exposure parameters
    sci_ngroups = output.meta.exposure.ngroups
    sci_nints = output.meta.exposure.nints

    # values defined for segmented data
    sci_int_start = output.meta.exposure.integration_start

    if sci_int_start is None:  # the data is not segmented
        sci_int_start = 1

    log.debug(f"RSCD correction using: nints={sci_nints}, ngroups={sci_ngroups}")
    log.debug(f"The first integration in the data is integration: {sci_int_start}")

    # For general RSCD flagging, we have to start with at least 4 groups. The last frame
    # has been rejected in the last frame correction, leaving us with 3 groups. We have to
    # have at least 2 valid groups to perform a fit. Therefore the minimum number of groups
    # to do an rscd flagging is 4 groups. MIRI has a set minimum of 5 groups (so only in rare
    # special cases will we have less than 4 groups).

    if sci_ngroups < 3:
        log.warning("Too few groups to apply RSCD correction")
        log.warning("RSCD step will be skipped")
        output.meta.cal_step.rscd = "SKIPPED"
        return output

    # Basic global checks:
    # ___________________
    # Will we have at least 3 groups left after skipping RSCD affected groups ? If not reduced
    # the number of groups skipped to ensure we have at least 3 groups. We need 3 because the
    # last frame correction has flagged the last group as DO_NOT_USE.

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

    # Note For segmented data the first integration in the file may not be the first
    # integration in the exposure. The value in meta.exposure.integration_start
    # holds the value of the first integration in the file.

    # Flag RSCD  groups in integration 1
    # __________________________________
    log.info(f"Number of groups to skip for integration 1: {group_skip_int1}")

    if sci_int_start == 1:  # we have segmented data and the data starts with the first integration.
        rscd_skip_array, num_rscd_lowered = flag_rscd(
            output,
            sci_int_start,
            sci_int_start,
            group_skip_int1,
            bright_use_2groups,
        )

        # print(rscd_skip_array.shape)

        output = apply_rscd_flags(output, sci_int_start, sci_int_start, rscd_skip_array)
        log.info(
            "Number of usable bright pixels with rscd flag groups "
            f"not set to DO_NOT_USE: {num_rscd_lowered}"
        )

        output.meta.rscd.keep_bright_firstgroup_int1 = num_rscd_lowered

        output.meta.rscd.ngroups_skip_int1 = group_skip_int1

    # Flag RSCD  groups in integration 2 and higher
    # ______________________________________________
    int_start = 1
    int_end = sci_nints - 1

    if sci_int_start != 1:  # we have segmented data and we are not on the first integration
        int_start = 0

    if sci_nints > 1:
        log.info(f"Number of groups to skip for integrations 2 and higher: {group_skip_int2p}")

        rscd_skip_array, num_rscd_lowered = flag_rscd(
            output, int_start, int_end, group_skip_int2p, bright_use_2groups
        )

        output = apply_rscd_flags(output, int_start, int_end, rscd_skip_array)
        log.info(
            "Number of usable bright pixels with rscd flag groups "
            f"not set to DO_NOT_USE: {num_rscd_lowered}"
        )

        output.meta.rscd.keep_bright_firstgroup_int2p = num_rscd_lowered
        output.meta.rscd.ngroups_skip_int2p = group_skip_int2p
    output.meta.cal_step.rscd = "COMPLETE"

    return output


def flag_rscd(output_model, int_start, int_end, rscd_skip, bright_use_2groups):
    """
    Find the initial groups to set to DO_NOT_USE based on RSCD rules.

    Parameters
    ----------
    output_model : RampModel
        Science data to be flagged

    int_start : int
        Starting integration

    int_end : int
        Endinging integration

    rscd_skip : int
        Number of groups to skip at the beginning of the ramp for integration range.

    bright_use_2groups : bool
        If True, for saturated data adjust number of rscd groups to flag to ensure we
        have at least 2 valid groups.

    Returns
    -------
    skip_array : np array int
        Array for each integration, group , pixel containing the number of rscd groups
        skip.
    num_rscd_lowered : int
        Number of pixels where the number of RSCD groups to flag as DO_NOT_USE was changed.
    """
    if int_end == 0:
        log.info(f"Number of groups to skip for integration 1 : {rscd_skip}")
    else:
        log.info(f"Number of groups to skip for integration 2 and higher : {rscd_skip}")

    n_ints = int_end - int_start + 1
    x_dim = output_model.groupdq.shape[3]
    y_dim = output_model.groupdq.shape[2]
    # print("x dim", x_dim)
    # print("y dim", y_dim)

    # print("number of ints", n_ints)
    skip_array = np.zeros((n_ints, y_dim, x_dim))
    skip_array[:, :, :] = rscd_skip

    if bright_use_2groups:
        # --- Logic for retaining at least 2 groups when we have saturating data
        min_group = rscd_skip + 2  # counting starting at 1

        is_sat_problem = (
            (
                output_model.groupdq[int_start : int_end + 1, min_group - 1, :, :]
                & dqflags.group["SATURATED"]
            )
            > 0
        ).astype(bool)

        # print(is_sat_problem.shape)
        num_sat = np.sum(is_sat_problem)
        num_rscd_lowered = num_sat

        log.info(
            f" There are {num_sat} saturated pixels that require the number of "
            "rscd groups flagged to be lowered"
        )
        if num_sat > 0:
            #  do dynamic rscd flagging - based on saturation group of every pixel
            x_dim = output_model.groupdq.shape[3]
            y_dim = output_model.groupdq.shape[2]
            # print("x dim", x_dim)
            # print("y dim", y_dim)

            n_ints = int_end - int_start + 1
            # print("number of ints", n_ints)
            skip_array = np.zeros((n_ints, y_dim, x_dim))
            skip_array[:, :, :] = rscd_skip
            # print("initial skip array shape", skip_array.shape)

            while num_sat > 0 and min_group > 0:
                skip_array[is_sat_problem] = np.maximum(skip_array[is_sat_problem] - 1, 0)
                # Collapse 3D to 2D: True if saturated in ANY integration
                is_sat_2d = np.any(is_sat_problem, axis=0)

                # Get Y and X coordinates
                y_coords, x_coords = np.where(is_sat_2d)

                output_model.pixeldq[y_coords, x_coords] |= dqflags.pixel["FLUX_ESTIMATED"]
                min_group = min_group - 1
                is_sat_problem = (
                    (
                        output_model.groupdq[int_start : int_end + 1, min_group - 1, :, :]
                        & dqflags.group["SATURATED"]
                    )
                    > 0
                ).astype(bool)

                num_sat = is_sat_problem.sum()

    return skip_array, num_rscd_lowered


def apply_rscd_flags(output_model, int_start, int_end, skip_array):
    """
    Find the initial groups to set to DO_NOT_USE based on RSCD rules.

    Parameters
    ----------
    output_model : RampModel
        Science data to be flagged

    int_start : int
        Starting integration

    int_end : int
        Endinging integration

    skip_array : int array
        Number of groups to skip at the beginning of the ramp for integration range.

    Returns
    -------
    output_model: RampModel
        Ramp datamodel with RSCD affected groups flagged as DO_NOT_USE
    """
    # 1. Extract the relevant region of the groupdq array
    # Shape: (N_ints, Groups, Y, X)
    dq = output_model.groupdq[int_start : int_end + 1, :, :, :]

    # 2. Create a grid of group indices
    # Shape: (Groups,) -> e.g., [0, 1, 2, 3...]
    num_groups = dq.shape[1]
    group_indices = np.arange(num_groups)

    # 3. Broadcast for comparison
    # We want: (1, Groups, 1, 1) < (N_ints, 1, Y, X)
    # This results in a 4D boolean mask
    mask = group_indices[None, :, None, None] < skip_array[:, None, :, :]

    # 4. Apply the DO_NOT_USE flag using the mask
    # This updates only the pixels/groups where the index is below the skip threshold
    dq[mask] |= dqflags.group["DO_NOT_USE"]

    # Put the modified dq back
    output_model.groupdq[int_start : int_end + 1, :, :, :] = dq

    return output_model


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
