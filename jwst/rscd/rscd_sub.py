#
#  Module for the RSCD correction for MIRI science data
#

import numpy as np
import logging
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, rscd_model, type):
    """
    Short Summary
    -------------
    The sole correction is to add the DO_NOT_USE flat to the GROUP data
    quality flags for the first N groups, where N is set by the rscd_model.

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd reference data

    type: string
        type of algorithm ['baseline' or 'enhanced']

    Returns
    -------
    output_model: data model object
        RSCD-corrected science data

    """

    # Save some data params for easy use later
    sci_ngroups = input_model.data.shape[1]

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # get the number of groups to flag
    nflag = get_rscd_parameters(input_model, rscd_model)

    # Update the step status, and if ngroups > nflag+3, set all of the GROUPDQ in
    # the first group to 'DO_NOT_USE'
    # the +3 is to ensure there is a slope to be fit including the flagging for
    # the last frame correction
    if sci_ngroups > (nflag + 3):
        output.groupdq[:, 0, :, :] = \
            np.bitwise_or(output.groupdq[:,0:nflag,:,:], dqflags.group['DO_NOT_USE'])
        log.debug(f"RSCD Sub: resetting GROUPDQ in the first {nflag} groups to DO_NOT_USE")
        output.meta.cal_step.rscd = 'COMPLETE'
    else:   # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.rscd = 'SKIPPED'

    return output


def get_rscd_parameters(input_model, rscd_model):

    """
    Read in the parameters from the reference file

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd reference data

    Returns
    -------
    nflag: int
        number of groups to flag

    """
    # determine the number of gropus to flag based on the subarray
    #   TBD

    # temp code
    nflag = 3

    return nflag
