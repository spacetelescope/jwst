import logging

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)

GOOD = dqflags.group["GOOD"]
DNU = dqflags.group["DO_NOT_USE"]
CHLO = dqflags.group["CHARGELOSS"]

CHLO_DNU = CHLO + DNU

__all__ = ["charge_migration"]


def charge_migration(output_model, signal_threshold):
    """
    Correct for charge migration.

    Each group in each ramp exceeding the signal threshold is flagged as
    CHARGELOSS and DO_NOT_USE, skipping groups already flagged as DO_NOT_USE.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        The input science data to be corrected. Updated in place.
    signal_threshold : float
        Science value above which a group will be flagged as CHARGELOSS.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        Data model with charge_migration applied; add CHARGELOSS and
        DO_NOT_USE flags to groups exceeding signal_threshold.
    """
    data = output_model.data
    gdq = output_model.groupdq

    log.info("Using signal_threshold: %.2f", signal_threshold)

    n_ints, n_grps, n_rows, n_cols = gdq.shape
    chargeloss_pix = np.where((data > signal_threshold) & (gdq != DNU))

    for k in range(len(chargeloss_pix[0])):
        integ, group = chargeloss_pix[0][k], chargeloss_pix[1][k]
        row, col = chargeloss_pix[2][k], chargeloss_pix[3][k]
        gdq[integ, group:, row, col] |= CHLO_DNU

        # North
        if row > 0:
            gdq[integ, group:, row - 1, col] |= CHLO_DNU

        # South
        if row < (n_rows - 1):
            gdq[integ, group:, row + 1, col] |= CHLO_DNU

        # East
        if col < (n_cols - 1):
            gdq[integ, group:, row, col + 1] |= CHLO_DNU

        # West
        if col > 0:
            gdq[integ, group:, row, col - 1] |= CHLO_DNU

    return output_model
