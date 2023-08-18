#  Module for undersampling correction
#
import logging
import numpy as np

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

GOOD = dqflags.group["GOOD"]
DNU = dqflags.group["DO_NOT_USE"]
UNSA = dqflags.group["UNDERSAMP"]
UNSA_DNU = UNSA + DNU


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
        Data model with undersampling_correction applied; add UNDERSAMP and
        DO_NOT_USE flags to groups exceeding signal_threshold
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
    Flag each group in each ramp that exceeds signal_threshold as UNDERSAMP and DO_NOT_USE,
    skipping groups already flagged as DO_NOT_USE.

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
    new_gdq : int, 4D array
        updated group dq array
    """
    n_ints, n_grps, n_rows, n_cols = gdq.shape
    ncols = data.shape[3]
    nrows = data.shape[2]

    new_gdq = gdq.copy()   # Updated gdq

    # Flag all exceedances with UNDERSAMP and NO_NOT_USE
    undersamp_pix = (data > signal_threshold) & (gdq != DNU)

    new_gdq[undersamp_pix] = np.bitwise_or(new_gdq[undersamp_pix], UNSA | DNU)

    # Reset groups previously flagged as DNU
    gdq_orig = gdq.copy()  # For resetting to previously flagged DNU
    wh_gdq_DNU = np.bitwise_and(gdq_orig, DNU)

    # Get indices for exceedances
    arg_where = np.argwhere(new_gdq == UNSA_DNU)
    a_int = arg_where[:, 0]  # array of integrations
    a_grp = arg_where[:, 1]  # array of groups
    a_row = arg_where[:, 2]  # array of rows
    a_col = arg_where[:, 3]  # array of columns

    # Process the 4 nearest neighbors of each exceedance
    # Pixel to the east
    xx_max_p1 = a_col[a_col < (ncols-1)] + 1
    i_int = a_int[a_col < (ncols-1)]
    i_grp = a_grp[a_col < (ncols-1)]
    i_row = a_row[a_col < (ncols-1)]

    if len(xx_max_p1) > 0:
        new_gdq[i_int, i_grp, i_row, xx_max_p1] = \
            np.bitwise_or(new_gdq[i_int, i_grp, i_row, xx_max_p1], UNSA | DNU)

    new_gdq[wh_gdq_DNU == 1] = gdq_orig[wh_gdq_DNU == 1]  # reset for earlier DNUs

    # Pixel to the west
    xx_m1 = a_col[a_col > 0] - 1
    i_int = a_int[a_col > 0]
    i_grp = a_grp[a_col > 0]
    i_row = a_row[a_col > 0]

    if len(xx_m1) > 0:
        new_gdq[i_int, i_grp, i_row, xx_m1] = \
            np.bitwise_or(new_gdq[i_int, i_grp, i_row, xx_m1], UNSA | DNU)

    new_gdq[wh_gdq_DNU == 1] = gdq_orig[wh_gdq_DNU == 1]  # reset for earlier DNUs

    # Pixel to the north
    yy_m1 = a_row[a_row > 0] - 1
    i_int = a_int[a_row > 0]
    i_grp = a_grp[a_row > 0]
    i_col = a_col[a_row > 0]

    if len(yy_m1) > 0:
        new_gdq[i_int, i_grp, yy_m1, i_col] = \
            np.bitwise_or(new_gdq[i_int, i_grp, yy_m1, i_col], UNSA | DNU)

    new_gdq[wh_gdq_DNU == 1] = gdq_orig[wh_gdq_DNU == 1]  # reset for earlier DNUs

    # Pixel to the south
    yy_max_p1 = a_row[a_row < (nrows-1)] + 1
    i_int = a_int[a_row < (nrows-1)]
    i_grp = a_grp[a_row < (nrows-1)]
    i_col = a_col[a_row < (nrows-1)]

    if len(yy_max_p1) > 0:
        new_gdq[i_int, i_grp, yy_max_p1, i_col] = \
            np.bitwise_or(new_gdq[i_int, i_grp, yy_max_p1, i_col], UNSA | DNU)

    new_gdq[wh_gdq_DNU == 1] = gdq_orig[wh_gdq_DNU == 1]  # reset for earlier DNUs

    return new_gdq
