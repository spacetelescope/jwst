
#  Module for undersampling correction
#
import logging
import numpy as np

from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def undersampling_correction(input_model, signal_threshold):
    """
    Short Summary
    -------------
    Correct for undersampling.

    Parameters
    ----------
    input_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected

    signal_threshold : float
        Science value above which a group will be flagged as DO_NOT_USE

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with undersampling_correction applied; add DO_NOT_USE flag
        to groups exceeding signal_threshold
    """
    data = input_model.data
    gdq = input_model.groupdq

    # Create the output model as a copy of the input
    output_model = input_model.copy()

    log.info('Using signal_threshold: %s', signal_threshold)

    gdq_new = flag_pixels(data, gdq, signal_threshold)

    # Save the flags in the output GROUPDQ array
    output_model.groupdq = gdq_new

    return output_model


def flag_pixels(data, gdq, signal_threshold):
    """
    Short Summary
    -------------
    Flag groups that exceed signal_threshold as DO_NOT_USE.

    Parameters
    ----------
    data : float, 4D array
        science array

    gdq : int, 4D array
        group dq array

    signal_threshold : float
        Science value above which a group will be flagged as DO_NOT_USE

    Returns
    -------
    gdq : int, 4D array
        updated group dq array
    """
    # For groups in the data that exceed the signal_threshold
    # add DO_NOT_USE to the group's DQ
    wh_exc = np.where(data >= signal_threshold)
    gdq[wh_exc] = np.bitwise_or(gdq[wh_exc], dqflags.group['DO_NOT_USE'])

    return gdq
