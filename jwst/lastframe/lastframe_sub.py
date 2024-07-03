#
#  Module for the lastframe correction for MIRI science data sets
#
import gc
import numpy as np
import logging

from stdatamodels.jwst.datamodels import dqflags

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

    # Create output as a copy of the input science data model
    output = input_model
    input_model.close()
    del input_model

    # Save some data params for easy use later
    sci_ngroups = output.data.shape[1]

    # Update the step status, and if ngroups > 2, set all of the GROUPDQ in
    # the final group to 'DO_NOT_USE'
    if sci_ngroups > 2:
        output.groupdq[:, -1, :, :] = \
            np.bitwise_or(output.groupdq[:, -1, :, :], dqflags.group['DO_NOT_USE'])
        log.debug("LastFrame Sub: resetting GROUPDQ in last frame to DO_NOT_USE")
        output.meta.cal_step.lastframe = 'COMPLETE'
    else:   # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.lastframe = 'SKIPPED'

    gc.collect()
    return output
