from __future__ import annotations
import warnings

import numpy as np

from jwst.datamodels import IFUImageModel  # type: ignore[attr-defined]
from stcal.outlier_detection.utils import medfilt
from stdatamodels.jwst.datamodels.dqflags import pixel


def badpix_selfcal(
    minimg: np.ndarray,
    flagfrac_lower: float = 0.001,
    flagfrac_upper: float = 0.001,
    kernel_size: int = 15,
    dispaxis=None,
) -> np.ndarray:
    """
    Flag residual artifacts as bad pixels in the DQ array of a JWST exposure.

    Parameters
    ----------
    minimg : np.ndarray
        Selfcal data of shape (x, y), i.e., after some operation has
        already been taken to combine multiple exposures,
        typically a MIN operation.
    flagfrac_lower : float
        Fraction of pixels to flag on the low end
    flagfrac_upper : float
        Fraction of pixels to flag on the high end
    kernel_size : int
        Size of kernel for median filter
    dispaxis : int
        Dispersion axis, either 1 or 2. If None, a two-dimensional
        median filter is applied.

    Returns
    -------
    flagged_indices : np.ndarray
        Indices of the flagged pixels,
        shaped like output from np.where
    """
    if dispaxis not in [1, 2, None]:
        raise ValueError("dispaxis must be either 1 or 2, or None.")

    # Determine outliers using median filter
    if dispaxis is None:
        kern_size = (kernel_size, kernel_size)
    elif dispaxis == 2:
        # check that this is the right way around!
        kern_size = (kernel_size, 1)
    elif dispaxis == 1:
        kern_size = (1, kernel_size)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")
        smoothed = medfilt(minimg, kern_size)
    minimg_hpf = minimg - smoothed

    # Flag outliers using percentile cutoff
    flag_low, flag_high = np.nanpercentile(
        minimg_hpf, [flagfrac_lower * 100, (1 - flagfrac_upper) * 100]
    )
    bad = (minimg_hpf > flag_high) | (minimg_hpf < flag_low)
    flagged_indices = np.where(bad)
    return flagged_indices


def apply_flags(input_model: IFUImageModel, flagged_indices: np.ndarray) -> IFUImageModel:
    """
    Apply the flagged indices to the input model.

    Sets the flagged pixels to NaN
    and the DQ flag to DO_NOT_USE + OTHER_BAD_PIXEL

    Parameters
    ----------
    input_model : IFUImageModel
        Input science data to be corrected
    flagged_indices : np.ndarray
        Indices of the flagged pixels,
        shaped like output from np.where

    Returns
    -------
    output_model : IFUImageModel
        Flagged data model
    """
    input_model.dq[flagged_indices] |= pixel["DO_NOT_USE"] + pixel["OTHER_BAD_PIXEL"]

    input_model.data[flagged_indices] = np.nan
    input_model.err[flagged_indices] = np.nan
    input_model.var_poisson[flagged_indices] = np.nan
    input_model.var_rnoise[flagged_indices] = np.nan
    input_model.var_flat[flagged_indices] = np.nan

    return input_model
