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
    Set data quality of the first group in an integration to DO_NOT_USE.

    This correction only works on MIRI data. If the number of groups is > 3,
    then the GROUP data quality flag of the first group is set to
    DO_NOT_USE.

    Parameters
    ----------
    output : DataModel
        Science data to be corrected
    bright_use_group1 : bool
        If True, setting first group data quality flag to DO_NOT_USE will not be
        done for pixels that have the saturation flag set for group 3.

    Returns
    -------
    output : DataModel
        Firstframe-corrected science data
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
                f"Number of bright pixels with first group not set to DO_NOT_USE, #{np.sum(svals)}"
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
