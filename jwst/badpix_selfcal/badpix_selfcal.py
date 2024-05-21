from __future__ import annotations
import numpy as np
import jwst.datamodels as dm
from jwst.outlier_detection.outlier_detection_ifu import medfilt


def badpix_selfcal(medbg: np.ndarray, 
                   flagfrac: float = 0.001, 
                   kernel_size: int = 15) -> np.ndarray:
    """
    Parameters
    ----------
    medbg: np.ndarray
        Background data of shape (x, y), i.e., after median has
        already been taken over the number of exposures
    flagfrac: float
        Fraction of pixels to flag on each of low and high end
    kernel_size: int
        Size of kernel for median filter

    Returns
    -------
    flagged_indices: np.ndarray
        Indices of the flagged pixels, 
        shaped like output from np.where
    """
    if flagfrac < 0 or flagfrac >=0.5:
        raise ValueError("flagfrac must be between 0 and 0.5. \
                        Note this fraction will be flagged on both high and low ends.")

    # Determine outliers using median filter
    kernel=np.array([1,kernel_size])
    smoothed = medfilt(medbg, kernel)
    medbg_hpf = medbg - smoothed

    # Flag outliers using percentile cutoff
    flag_low, flag_high = np.nanpercentile(medbg_hpf, [flagfrac*100, (1-flagfrac)*100])
    bad = (medbg_hpf > flag_high) | (medbg_hpf < flag_low)
    flagged_indices = np.where(bad)

    return flagged_indices


def apply_flags(input_model: dm.IFUImageModel, flagged_indices: np.ndarray) -> dm.IFUImageModel:
    """
    Apply the flagged indices to the input model.

    Parameters
    ----------
    input_model: dm.IFUImageModel
        Input science data to be corrected
    flagged_indices: np.ndarray
        Indices of the flagged pixels, 
        shaped like output from np.where

    Returns
    -------
    output_model: dm.IFUImageModel
        Data model with flagged pixels set to NaN
        in data and err arrays, and dq flag set to 1
    """

    input_model.data[flagged_indices] = np.nan
    input_model.err[flagged_indices] = np.nan
    input_model.dq[flagged_indices] = 1

    return input_model
