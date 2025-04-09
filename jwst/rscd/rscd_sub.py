#
#  Module for the RSCD correction for MIRI science data
#

import numpy as np
import logging

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(output_model, rscd_model):
    """
    Set the initial groups of MIRI data to 'DO_NOT_USE' for integrations 2 and higher.

    The number of initial groups to set to 'DO_NOT_USE' is read in from the RSCD reference
    file. No flagging is done for the first integration.

    Parameters
    ----------
    output_model : RampModel
        Input ramp datamodel

    rscd_model : RSCDModel
        RSCD reference datamodel

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

    group_skip = param["skip"]
    output_model = correction_skip_groups(output_model, group_skip)
    return output_model


def correction_skip_groups(output, group_skip):
    """
    Set the initial groups in integration to DO_NOT_USE to skip groups affected by RSCD effect.

    Parameters
    ----------
    output : RampModel
        Science data to be flagged

    group_skip : int
        Number of groups to skip at the beginning of the ramp

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
    log.info(f"Number of groups to skip for integrations 2 and higher: {group_skip}")

    # If ngroups <= group_skip+3, skip the flagging
    # the +3 is to ensure there is a slope to be fit including the flagging for
    # the last frame correction
    if sci_ngroups <= (group_skip + 3):
        log.warning("Too few groups to apply RSCD correction")
        log.warning("RSCD step will be skipped")
        output.meta.cal_step.rscd = "SKIPPED"
        return output

    # The RSCD correction is applied to integrations 2 and higher.
    # For segmented data the first integration in the file may not be the first
    # integration in the exposure. The value in meta.exposure.integration_start
    # holds the value of the first integration in the file.
    # If a correction is to be done and if ngroups > group_skip+3, then set all
    # of the GROUPDQ flags from groups 0 to group_skip to 'DO_NOT_USE'.

    int_start = 1
    if sci_int_start != 1:  # we have segmented data
        int_start = 0

    output.groupdq[int_start:, 0:group_skip, :, :] = np.bitwise_or(
        output.groupdq[int_start:, 0:group_skip, :, :], dqflags.group["DO_NOT_USE"]
    )
    log.debug(f"RSCD Sub: adding DO_NOT_USE to GROUPDQ for the first {group_skip} groups")
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
        group_skip_table = tabdata["group_skip"]
        if subarray_table == subarray and readpatt_table == readpatt:
            param["skip"] = group_skip_table
            break

    return param
