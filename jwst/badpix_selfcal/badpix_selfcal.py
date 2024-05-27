from __future__ import annotations
import numpy as np
import jwst.datamodels as dm
from jwst.outlier_detection.outlier_detection_ifu import medfilt
from stdatamodels.jwst.datamodels.dqflags import pixel


def badpix_selfcal(medbg: np.ndarray, 
                   flagfrac: float = 0.001, 
                   kernel_size: int = 15) -> np.ndarray:
    """
    Flag residual artifacts as bad pixels in the DQ array of a JWST exposure

    Parameters
    ----------
    medbg: np.ndarray
        Background data of shape (x, y), i.e., after some operation has
        already been taken to combine multiple exposures,
        typically a MIN operation.
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
    Apply the flagged indices to the input model. Sets the flagged pixels to NaN
    and the DQ flag to DO_NOT_USE + OTHER_BAD_PIXEL

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
        Flagged data model
    """

    input_model.dq[flagged_indices] = pixel["DO_NOT_USE"] + pixel["OTHER_BAD_PIXEL"]

    input_model.data[flagged_indices] = np.nan
    input_model.err[flagged_indices] = np.nan
    input_model.var_poisson[flagged_indices] = np.nan
    input_model.var_rnoise[flagged_indices] = np.nan
    input_model.var_flat[flagged_indices] = np.nan

    return input_model
