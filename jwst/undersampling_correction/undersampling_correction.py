#  Module for undersampling correction
#
import logging
import numpy as np

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def undersampling_correction(input_model, signal_threshold):
    """
    Correct for undersampling.

    Parameters
    ----------
    input_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected

    signal_threshold : float
        Science value above which a group will be flagged as UNDERSAMP

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with undersampling_correction applied; add UNDERSAMP flag
        to groups exceeding signal_threshold
    """
    data = input_model.data
    gdq = input_model.groupdq

    # Create the output model as a copy of the input
    output_model = input_model.copy()

    log.info('Using signal_threshold: %.2f', signal_threshold)

    gdq_new = flag_pixels(data, gdq, signal_threshold)

    # Save the flags in the output GROUPDQ array
    output_model.groupdq = gdq_new

    return output_model


def flag_pixels(data, gdq, signal_threshold):
    """
    Flag first group in each ramp that exceeds signal_threshold as UNDERSAMP and DO_NOT_USE,
    skipping groups already flagged as DO_NOT_USE; then flag all subsequent groups in the ramp.

    Parameters
    ----------
    data : float, 4D array
        science array

    gdq : int, 4D array
        group dq array

    signal_threshold : float
        Science value above which a group will be flagged as UNDERSAMP and DO_NOT_USE

    Returns
    -------
    gdq : int, 4D array
        updated group dq array
    """
    n_ints, n_grps, n_rows, n_cols = gdq.shape
    num_pix = n_cols * n_rows

    lowest_exc_1d = np.zeros(num_pix) + n_grps

    for ii_int in range(n_ints):
        for ii_grp in range(n_grps):
            data_1d = data[ii_int, ii_grp, :, :].reshape(num_pix)  # vectorize slice
            gdq_1d = gdq[ii_int, ii_grp, :, :].reshape(num_pix)

            wh_not_dnu = np.logical_not(gdq_1d & dqflags.group['DO_NOT_USE'])

            # In the current group for all ramps, locate pixels that :
            #  a) exceed the signal_threshold, and
            #  b) have not been previously flagged as an exceedance, and
            #  c) were not flagged in an earlier step as DO_NOT_USE
            wh_exc_1d = np.where((data_1d > signal_threshold) &
                                 (lowest_exc_1d == n_grps) & wh_not_dnu)

            # ... and mark those pixels, as current group is their first exceedance
            if len(wh_exc_1d[0] > 0):  # For ramps previously unflagged ...
                lowest_exc_1d[wh_exc_1d] = ii_grp

    # Flag current and subsequent groups
    lowest_exc_2d = lowest_exc_1d.reshape((n_rows, n_cols))
    for ii_int in range(n_ints):
        for ii_grp in range(n_grps):
            wh_set_flag = np.where(lowest_exc_2d == ii_grp)

            # set arrays of components
            yy = wh_set_flag[0]
            xx = wh_set_flag[1]

            gdq[ii_int, ii_grp:, yy, xx] = \
                np.bitwise_or(gdq[ii_int, ii_grp:, yy, xx], dqflags.group['UNDERSAMP']
                              | dqflags.group['DO_NOT_USE'])

    return gdq
