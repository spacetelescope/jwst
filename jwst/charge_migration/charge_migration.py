#  Module for charge migration
#
import logging
import numpy as np

from stdatamodels.jwst.datamodels import dqflags


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

GOOD = dqflags.group["GOOD"]
DNU = dqflags.group["DO_NOT_USE"]
CHLO = dqflags.group["CHARGELOSS"]

CHLO_DNU = CHLO + DNU


def charge_migration(output_model, signal_threshold):
    """
    Correct for charge migration.

    Parameters
    ----------
    output_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected.
    signal_threshold : float
        Science value above which a group will be flagged as CHARGELOSS.

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with charge_migration applied; add CHARGELOSS and
        DO_NOT_USE flags to groups exceeding signal_threshold.
    """
    data = output_model.data
    gdq = output_model.groupdq

    log.info("Using signal_threshold: %.2f", signal_threshold)

    gdq_new = flag_pixels(data, gdq, signal_threshold)

    # Save the flags in the output GROUPDQ array
    output_model.groupdq = gdq_new

    return output_model


def flag_pixels(data, gdq, signal_threshold):
    """
    Flag each group in each ramp that exceeds signal_threshold.

    Each group in each ramp exceeding the signal threshold is flagged as
    CHARGELOSS and DO_NOT_USE, skipping groups already flagged as DO_NOT_USE.

    Parameters
    ----------
    data : float, 4D array
        Science array.
    gdq : int, 4D array
        Group DQ array.
    signal_threshold : float
        Science value above which a group will be flagged as CHARGELOSS and DO_NOT_USE.

    Returns
    -------
    new_gdq : int, 4D array
        Updated group DQ array.
    """
    n_ints, n_grps, n_rows, n_cols = gdq.shape
    chargeloss_pix = np.where((data > signal_threshold) & (gdq != DNU))

    new_gdq = gdq.copy()

    for k in range(len(chargeloss_pix[0])):
        integ, group = chargeloss_pix[0][k], chargeloss_pix[1][k]
        row, col = chargeloss_pix[2][k], chargeloss_pix[3][k]
        new_gdq[integ, group:, row, col] |= CHLO_DNU

        # North
        if row > 0:
            new_gdq[integ, group:, row - 1, col] |= CHLO_DNU

        # South
        if row < (n_rows - 1):
            new_gdq[integ, group:, row + 1, col] |= CHLO_DNU

        # East
        if col < (n_cols - 1):
            new_gdq[integ, group:, row, col + 1] |= CHLO_DNU

        # West
        if col > 0:
            new_gdq[integ, group:, row, col - 1] |= CHLO_DNU

    return new_gdq
