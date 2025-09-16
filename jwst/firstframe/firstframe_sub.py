#
#  Module for the firstframe correction for MIRI science data sets
#
import logging

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)

__all__ = ["do_correction"]


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
            tvals[~svals] |= dqflags.group["DO_NOT_USE"]

            # also set a FLUX_ESTIMATED flag in the 2D pixeldq image
            # for any group1 pixels that were kept
            any_kept = np.any(svals, axis=0)
            output.pixeldq[any_kept] |= dqflags.pixel["FLUX_ESTIMATED"]

            any_good_kept = (output.pixeldq[any_kept] & dqflags.pixel["DO_NOT_USE"]) == 0
            log.info(
                "Number of usable bright pixels with first group "
                f"not set to DO_NOT_USE: {np.sum(any_good_kept)}"
            )
        else:
            output.groupdq[:, 0, :, :] |= dqflags.group["DO_NOT_USE"]

        log.debug("FirstFrame Sub: resetting GROUPDQ in first frame to DO_NOT_USE")
        output.meta.cal_step.firstframe = "COMPLETE"
    else:  # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.firstframe = "SKIPPED"

    return output
