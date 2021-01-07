#
#  Module for the lastframe correction for MIRI science data sets
#
import numpy as np
import logging
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model):
    """
    Short Summary
    -------------
    The sole correction is to reset to DO_NOT_USE the GROUP data quality flags
    for the final group, if the number of groups is greater than 2.

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    Returns
    -------
    output: data model object
        lastframe-corrected science data

    """

    # Save some data params for easy use later
    sci_ngroups = input_model.data.shape[1]

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # Update the step status, and if ngroups > 2, set all of the GROUPDQ in
    # the final group to 'DO_NOT_USE'
    if sci_ngroups > 2:
        output.groupdq[:,-1,:,:] = \
            np.bitwise_or(output.groupdq[:,-1,:,:], dqflags.group['DO_NOT_USE'])
        log.debug("LastFrame Sub: resetting GROUPDQ in last frame to DO_NOT_USE")
        output.meta.cal_step.lastframe = 'COMPLETE'
    else:   # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.lastframe = 'SKIPPED'

    return output
