import numpy as np
from jwst.resample.resample_utils import build_mask

from jwst import datamodels as dm

from stcal.outlier_detection.utils import compute_weight_threshold
from .utils import flag_model_crs, nanmedian3D
from ._fileio import save_median

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["detect_outliers"]


def detect_outliers(
    input_model,
    save_intermediate_results,
    good_bits,
    maskpt,
    rolling_window_width,
    snr,
    make_output_path,
):
    """
    Flag outliers in tso data.

    Parameters
    ----------
    input_model : ~jwst.datamodels.CubeModel
        The input cube model.
    save_intermediate_results : bool
        If True, save the rolling median model as a CubeModel.
    good_bits : int
        DQ flag bit values indicating good pixels.
    maskpt : float
        The percentage of the mean weight to use as a threshold for masking.
    rolling_window_width : int
        The width of the rolling median window.
    snr : float
        The signal-to-noise ratio threshold for flagging outliers.
    make_output_path : callable
        A function that generates a path for saving intermediate results.

    Returns
    -------
    ~jwst.datamodels.CubeModel
        The input model with outliers flagged.
    """
    if not isinstance(input_model, dm.JwstDataModel):
        input_model = dm.open(input_model)
    if isinstance(input_model, dm.ModelContainer):
        raise TypeError("OutlierDetectionTSO does not support ModelContainer input.")
    weighted_cube = weight_no_resample(input_model, good_bits)

    weight_threshold = compute_weight_threshold(weighted_cube.wht, maskpt)

    if (rolling_window_width > 1) and (rolling_window_width < weighted_cube.shape[0]):
        medians = compute_rolling_median(weighted_cube, weight_threshold, w=rolling_window_width)

    else:
        medians = nanmedian3D(weighted_cube.data, overwrite_input=False)
        # this is a 2-D array, need to repeat it into the time axis
        # for consistent shape with rolling median case
        medians = np.broadcast_to(medians, weighted_cube.shape)

    # Save median model if pars['save_intermediate_results'] is True
    # this will be a CubeModel with rolling median values.
    if save_intermediate_results:
        median_model = dm.CubeModel(data=medians)  # type: ignore[name-defined]
        with dm.open(weighted_cube) as dm0:
            median_model.update(dm0)
        save_median(median_model, make_output_path)
        del median_model

    # no need for blotting, resample is turned off for TSO
    # go straight to outlier detection
    log.info("Flagging outliers")
    flag_model_crs(
        input_model,
        medians,
        snr,
    )
    return input_model


def weight_no_resample(input_model, good_bits):
    """
    Give weights to model without resampling.

    Parameters
    ----------
    input_model : ~jwst.datamodels.CubeModel
        The input cube model.
    good_bits : int
        DQ flag bit values indicating good pixels.

    Returns
    -------
    ~jwst.datamodels.CubeModel
        A copy of the input cube model with weights assigned in the `wht` extension.

    Notes
    -----
    Prior to PR #8473, the `build_driz_weight` function was used to
    create the weights for the input models for TSO data. However, that
    function was simply returning a copy of the DQ array because the
    var_noise was not being passed in by calwebb_tso3. As of PR #8473,
    a cube model that includes the var_noise is passed into TSO
    outlier detection, so `build_driz_weight` would weight the cube model
    by the variance. Therefore `build_driz_weight` was removed in order to
    preserve the original behavior. If it is determined later that exposure
    time or inverse variance weighting should be used here, build_driz_weight
    should be re-implemented.
    """
    weighted_cube = input_model.copy()
    dqmask = build_mask(input_model.dq, good_bits)
    weighted_cube.wht = dqmask.astype(np.float32)
    return weighted_cube


def compute_rolling_median(
    model: dm.CubeModel,  # type: ignore[name-defined]
    weight_threshold: np.ndarray,
    w: int = 25,
) -> np.ndarray:
    """
    Set bad and low-weight data to NaN, then compute the rolling median over the time axis.

    Parameters
    ----------
    model : ~jwst.datamodels.CubeModel
        The input cube model

    weight_threshold : np.ndarray
        The weight thresholds for each integration.

    w : int
        The window size for the rolling median.

    Returns
    -------
    np.ndarray
        The rolling median of the input data. Same dimensions as input.
    """
    sci = model.data
    weight = model.wht
    badmask = np.less(weight, weight_threshold)
    log.debug(f"Percentage of pixels with low weight: {np.sum(badmask) / weight.size * 100}")

    # Fill resampled_sci array with nan's where mask values are True
    sci[badmask] = np.nan
    del badmask

    if w > sci.shape[0]:
        raise ValueError("Window size must be less than the number of integrations.")
    meds = moving_median_over_zeroth_axis(sci, w)

    del sci
    return meds


def moving_median_over_zeroth_axis(x: np.ndarray, w: int) -> np.ndarray:
    """
    Calculate the median of a moving window over the zeroth axis of an N-d array.

    Algorithm works by expanding the array into an additional dimension
    where the new axis has the same length as the window size. Each entry in that
    axis is a copy of the original array shifted by 1 with respect to the previous
    entry, such that the rolling median is simply the median over the new axis.
    modified from https://stackoverflow.com/a/71154394, see link for more details.

    Parameters
    ----------
    x : np.ndarray
        The input array.
    w : int
        The window size.

    Returns
    -------
    np.ndarray
        The rolling median of the input array. Same dimensions as input.
    """
    if w <= 1:
        raise ValueError("Rolling median window size must be greater than 1.")
    shifted = np.zeros((x.shape[0] + w - 1, w, *x.shape[1:])) * np.nan
    for idx in range(w - 1):
        shifted[idx : -w + idx + 1, idx] = x
    shifted[idx + 1 :, idx + 1] = x
    medians: np.ndarray = np.median(shifted, axis=1)
    for idx in range(w - 1):
        medians[idx] = np.median(shifted[idx, : idx + 1])
        medians[-idx - 1] = np.median(shifted[-idx - 1, -idx - 1 :])
    medians = medians[(w - 1) // 2 : -(w - 1) // 2]

    # Fill in the edges with the nearest valid value
    medians[: w // 2] = medians[w // 2]
    medians[-w // 2 :] = medians[-w // 2]
    return medians
