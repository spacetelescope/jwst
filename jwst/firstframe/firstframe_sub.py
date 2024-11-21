#
#  Module for the firstframe correction for MIRI science data sets
#
import numpy as np
import logging

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(output, bright_use_group1=False):
    """
    Short Summary
    -------------
    The sole correction is to add the DO_NOT_USE flat to the GROUP data
    quality flags for the first group, if the number of groups is greater
    than 3.

    Parameters
    ----------
    output: data model object
        science data to be corrected
    bright_use_group1: boolean
        do not flag group1 for bright pixels = group 3 saturated

    Returns
    -------
    output: data model object
        firstframe-corrected science data

    """

    # Save some data params for easy use later
    sci_ngroups = output.data.shape[1]

    # Update the step status, and if ngroups > 3, set all GROUPDQ in
    # the first group to 'DO_NOT_USE'
    if sci_ngroups > 3:
        if bright_use_group1:
            # do not set DO_NOT_USE in the case where saturation happens in
            # group3 as in this case the first frame effect is small compared to the
            # signal in group2-group1
            svals = (output.groupdq[:, 2, :, :] & dqflags.group["SATURATED"]) > 0
            tvals = output.groupdq[:, 0, :, :]
            tvals[~svals] = np.bitwise_or(
                (output.groupdq[:, 0, :, :])[~svals], dqflags.group["DO_NOT_USE"]
            )
            output.groupdq[:, 0, :, :] = tvals
            log.info(
                f"FirstFrame Sub: bright_first_frame set, #{np.sum(svals)} bright pixels group1 not set to DO_NOT_USE"
            )
        else:
            output.groupdq[:, 0, :, :] = np.bitwise_or(
                output.groupdq[:, 0, :, :], dqflags.group["DO_NOT_USE"]
            )

        log.debug("FirstFrame Sub: resetting GROUPDQ in first frame to DO_NOT_USE")
        output.meta.cal_step.firstframe = "COMPLETE"
    else:  # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.firstframe = "SKIPPED"

    return output
