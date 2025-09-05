import logging
import warnings
from fractions import Fraction

import numpy as np
from astropy.stats import sigma_clipped_stats
from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise import clean_flicker_noise as cfn
from jwst.lib.basic_utils import disable_logging

log = logging.getLogger(__name__)

__all__ = ["quick_clean", "niriss_image_parameters", "nircam_image_parameters"]


def quick_clean(input_model, flat_filename=None):
    """
    Run a quick version of the flicker noise cleaning on a rate image.

    This will return a median-cleaned and optionally flat-fielded rate image.
    Input data may be a ramp, rate, or rateints model.  If rateints, the
    data from all integrations are averaged before cleaning.

    All cleaning parameters are set to default values.

    Parameters
    ----------
    input_model : RampModel, ImageModel, or CubeModel
        Input data to clean.
    flat_filename : str or None
        If provided, will be read in and divided into the image before cleaning.

    Returns
    -------
    image : ndarray of float
        2D cleaned rate image.
    mask : ndarray of bool
        2D mask marking background pixels as True.
    """
    # Make a draft rate file first
    if isinstance(input_model, datamodels.RampModel):
        rate_model = cfn.make_rate(input_model)
    elif isinstance(input_model, datamodels.CubeModel):
        # Input is a rateints file:
        # just average the data across the integrations
        rate_model = datamodels.ImageModel(np.nanmean(input_model.data, axis=0))
        rate_model.update(input_model)
    else:
        # Input is already a rate file
        rate_model = datamodels.ImageModel(input_model.data.copy())
        rate_model.update(input_model)

    # Image to clean
    image = rate_model.data

    # Divide by flat
    flat = cfn.read_flat_file(rate_model, flat_filename)
    if flat is not None:
        image /= flat

    # Make a scene mask
    mask = np.full(image.shape, True)
    cfn.clip_to_background(image, mask)

    # Find and replace/mask NaNs
    nan_pix = np.isnan(image)
    image[nan_pix] = 0.0
    mask[nan_pix] = False

    # Subtract median background
    background = cfn.background_level(image, mask)
    image -= background

    # Flag more signal in the background subtracted image,
    # with sigma set by the lower half of the distribution only
    cfn.clip_to_background(image, mask, lower_half_only=True)

    # Median-clean the background subtracted image
    axis_to_correct = abs(rate_model.meta.subarray.slowaxis)
    image = cfn.median_clean(image, mask, axis_to_correct)

    # Restore NaNs and background to cleaned image
    image += background
    image[nan_pix] = np.nan

    rate_model.close()
    del rate_model
    return image, mask


def _row_average(image, start_fraction=0.0, stop_fraction=1.0):
    """
    Average an image section over specified rows for all columns.

    Parameters
    ----------
    image : ndarray
        Full 2D image.
    start_fraction : float or Fraction, optional
        Start index for the rows to average over.
    stop_fraction : float or Fraction, optional
        Stop index for the rows to average over.

    Returns
    -------
    ndarray
        1D array with size image.shape[0], containing average values for each column.
    """
    shape = image.shape
    y_start = int(round(start_fraction * shape[0]))
    y_stop = int(round(stop_fraction * shape[0]))
    return np.nanmean(image[y_start:y_stop, :], axis=0)


def _col_average(image, start_fraction=0.0, stop_fraction=1.0):
    """
    Average an image section over specified columns for all rows.

    Parameters
    ----------
    image : ndarray
        Full 2D image.
    start_fraction : float or Fraction, optional
        Start index for the columns to average over.
    stop_fraction : float or Fraction, optional
        Stop index for the columns to average over.

    Returns
    -------
    ndarray
        1D array with size image.shape[1], containing average values for each column.
    """
    shape = image.shape
    x_start = int(round(start_fraction * shape[1]))
    x_stop = int(round(stop_fraction * shape[1]))
    return np.nanmean(image[:, x_start:x_stop], axis=1)


def _line_fit(data):
    """
    Fit a line to 1D data.

    Parameters
    ----------
    data : ndarray
        1D data to fit.

    Returns
    -------
    intercept : float
        Constant coefficient for the fit line.
    slope : float
        Slope coefficient for the fit line.
    sigma : float
        Standard deviation of the residuals on the fit.
    """
    good = np.isfinite(data)
    x = np.arange(data.size)
    p_fit = np.polynomial.Polynomial.fit(x[good], data[good], deg=1)

    # Get the line coefficients
    coef = p_fit.convert().coef
    intercept = coef[0]
    if len(coef) > 1:
        slope = coef[1]
    else:
        slope = 0.0

    # Get the standard deviation of the residuals on the fit
    _, _, sigma = sigma_clipped_stats(data[good] - p_fit(x[good]), sigma=10.0)

    return intercept, slope, sigma


def _max_channel_offset(data):
    """
    Calculate the maximum average offset between channels.

    If the data is not full frame, returns 0.0.

    Parameters
    ----------
    data : ndarray
        1D array of average values, taken across the channels.
        If size is not 2048, the maximum offset is assumed to be 0.0.

    Returns
    -------
    max_offset : float
        The absolute maximum value in
        [channel2 - channel1, channel3 - channel2, channel4 - channel 3].
    """
    if data.size != 2048:
        return 0.0

    n_output = 4
    channel_size = 512
    average_values = np.zeros(n_output)
    cstart = 0
    cstop = channel_size
    for channel in range(n_output):
        av, _, _ = sigma_clipped_stats(data[cstart:cstop], sigma=3.0)
        average_values[channel] = av
        cstart += channel_size
        cstop += channel_size

    max_offset = np.max(np.abs(average_values[1:] - average_values[:-1]))
    return max_offset


def _image_decision_tree(
    parameters, limits, mask_frac, max_offset, max_slope, slope_ratio, channel_sigma
):
    """
    Decide parameter values from input statistics.

    Parameters are updated in place.  Only "fit_by_channel" and "background_method"
    will be modified.

    Parameters
    ----------
    parameters : dict
        Parameters to override.
    limits : dict
        Heuristic limits for statistics in various combinations.
    mask_frac : float
        Mask fraction value to check against limits.
    max_offset : float
        Max channel offset value to check against limits.
    max_slope : float
        Max slope value to check against limits.
    slope_ratio : float
        Slope ratio value to check against limits.
    channel_sigma : float
        Channel standard deviation value to check against limits.
    """
    log.debug("Image stats:")
    log.debug(f"  mask fraction: {mask_frac:.5g}")
    log.debug(f"  max offset: {max_offset:.5g}")
    log.debug(f"  max slope: {max_slope:.5g}")
    log.debug(f"  slope ratio: {slope_ratio:.5g}")
    log.debug(f"  channel sigma: {channel_sigma:.5g}")
    log.debug("Autoparam decision tree:")

    # Check input against limits
    high_mask_frac = mask_frac >= limits["mask_frac"]
    high_max_offset = max_offset >= limits["max_offset_high_mask"]
    high_max_slope_low_offset = max_slope >= limits["max_slope_low_offset"]
    high_max_slope_high_offset = max_slope >= limits["max_slope_high_offset"]
    high_slope_ratio = slope_ratio >= limits["slope_ratio"]
    high_channel_sigma = channel_sigma >= limits["channel_sigma"]

    if not high_mask_frac and high_max_offset:
        log.debug("  Low mask, high offset: median background, fit by channel")
        parameters["background_method"] = "median"
        parameters["fit_by_channel"] = True
    elif high_mask_frac:
        parameters["fit_by_channel"] = False
        if not high_max_offset:
            if high_max_slope_low_offset:
                log.debug("  High mask, low offset, high slope: model background, fit whole image")
                parameters["background_method"] = "model"
            else:
                log.debug("  High mask, low offset, low slope: median background, fit whole image")
                parameters["background_method"] = "median"
        else:
            if high_max_slope_high_offset:
                log.debug("  High mask, high offset, high slope: model background, fit whole image")
                parameters["background_method"] = "model"
            else:
                log.debug("  High mask, high offset, low slope: median background, fit whole image")
                parameters["background_method"] = "median"
    else:
        if high_max_slope_high_offset:
            log.debug("  Low mask, low offset, high slope: model background")
            parameters["background_method"] = "model"
        elif high_max_slope_low_offset:
            if not high_slope_ratio:
                log.debug("  Low mask, low offset, medium slope, low slope ratio: model background")
                parameters["background_method"] = "model"
            else:
                log.debug(
                    "  Low mask, low offset, medium slope, high slope ratio: median background"
                )
                parameters["background_method"] = "median"
        else:
            log.debug("  Low mask, low offset, low slope: median background")
            parameters["background_method"] = "median"
        if max_offset > 0 and high_channel_sigma:
            log.debug("  High channel sigma: fit by channel")
            parameters["fit_by_channel"] = True
        else:
            log.debug("  Low channel sigma: fit whole image")
            parameters["fit_by_channel"] = False


def niriss_image_parameters(input_model, flat_filename):
    """
    Determine appropriate parameters for cleaning a NIRISS image.

    Parameters
    ----------
    input_model : RampModel, ImageModel, or CubeModel
        Input model.
    flat_filename : str
        Path to a FLAT reference file.

    Returns
    -------
    dict
        Step parameters to override.
    """
    # Check for a flat file
    if flat_filename is None or flat_filename == "N/A":
        apply_flat = False
        flat_filename = None
    else:
        apply_flat = True

    # Set up some baseline defaults
    parameters = {
        "apply_flat_field": apply_flat,
        "background_method": "median",
        "fit_by_channel": False,
    }

    # Make a rate file and run a quick default median clean on it.
    with disable_logging(logging.ERROR):
        cleaned_image, mask = quick_clean(input_model, flat_filename)

    # Fit a 2D background to the cleaned image
    background_2d = cfn.background_level(cleaned_image, mask, background_method="model")

    # Replace sources with background values
    if np.isscalar(background_2d):
        cleaned_image[~mask] = background_2d
    else:
        cleaned_image[~mask] = background_2d[~mask]

    # Get the masked fraction
    mask_frac = np.sum(~mask) / mask.size

    # Get average row and column values from the middle quarter
    # of the cleaned background image
    avg_over_row = _row_average(
        cleaned_image, start_fraction=Fraction(3, 8), stop_fraction=Fraction(5, 8)
    )
    avg_over_col = _col_average(
        cleaned_image, start_fraction=Fraction(3, 8), stop_fraction=Fraction(5, 8)
    )

    # Fit a line to the column and row stats to get slopes at the
    # center of the image
    _, col_slope, _ = _line_fit(avg_over_row)
    _, row_slope, _ = _line_fit(avg_over_col)

    # Get the maximum slope and ratio of min to max
    abs_slopes = np.abs([row_slope, col_slope])
    max_slope = np.max(abs_slopes)
    min_slope = np.min(abs_slopes)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "divide by zero", RuntimeWarning)
        slope_ratio = max_slope / min_slope

    # Also get an average from rows near the top of the image:
    # for NIRISS, channel 1 can look significantly different from the others
    ch1_avg_over_row = _row_average(
        cleaned_image, start_fraction=Fraction(10, 12), stop_fraction=Fraction(11, 12)
    )

    # Get the standard deviation of the residuals on the fit for channel 1
    _, _, ch1_sigma = _line_fit(ch1_avg_over_row)

    # From the average over column levels, get the maximum channel offset
    max_offset = _max_channel_offset(avg_over_col)

    # Use a heuristic decision tree from calculated stats
    # to update the override parameters in place
    limits = {
        "mask_frac": 0.07,
        "max_offset_low_mask": 0.02,
        "max_offset_high_mask": 0.07,
        "max_slope_low_offset": 1e-5,
        "max_slope_high_offset": 1e-4,
        "slope_ratio": 5.0,
        "channel_sigma": 0.02,
    }
    _image_decision_tree(
        parameters, limits, mask_frac, max_offset, max_slope, slope_ratio, ch1_sigma
    )

    return parameters


def nircam_image_parameters(input_model, flat_filename):
    """
    Determine appropriate parameters for cleaning a NIRCam image.

    Parameters
    ----------
    input_model : RampModel, ImageModel, or CubeModel
        Input model.
    flat_filename : str
        Path to a FLAT reference file.

    Returns
    -------
    dict
        Step parameters to override.
    """
    # Check for a flat file
    if flat_filename is None or flat_filename == "N/A":
        apply_flat = False
        flat_filename = None
    else:
        apply_flat = True

    # Set up some baseline defaults
    parameters = {
        "apply_flat_field": apply_flat,
        "background_method": "median",
        "fit_by_channel": False,
    }

    # Make a rate file and run a quick default median clean on it.
    with disable_logging(logging.ERROR):
        cleaned_image, mask = quick_clean(input_model, flat_filename)

    # Fit a 2D background to the cleaned image
    background_2d = cfn.background_level(cleaned_image, mask, background_method="model")

    # Replace sources with background values
    if np.isscalar(background_2d):
        cleaned_image[~mask] = background_2d
    else:
        cleaned_image[~mask] = background_2d[~mask]

    # Get the masked fraction
    mask_frac = np.sum(~mask) / mask.size

    # Get average row and column values from the middle quarter
    # of the cleaned image
    avg_over_row = _row_average(
        cleaned_image, start_fraction=Fraction(3, 8), stop_fraction=Fraction(5, 8)
    )
    avg_over_col = _col_average(
        cleaned_image, start_fraction=Fraction(3, 8), stop_fraction=Fraction(5, 8)
    )

    # Fit a line to the column and row stats to get slopes at the center.
    # From the column fit, also get the standard deviation of the
    # fit residuals across channels.
    _, col_slope, channel_sigma = _line_fit(avg_over_row)
    _, row_slope, _ = _line_fit(avg_over_col)

    # Get the maximum slope and ratio of min to max
    abs_slopes = np.abs([row_slope, col_slope])
    max_slope = np.max(abs_slopes)
    min_slope = np.min(abs_slopes)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "divide by zero", RuntimeWarning)
        slope_ratio = max_slope / min_slope

    # From the average column levels, get the maximum channel offset
    max_offset = _max_channel_offset(avg_over_row)

    # Use a heuristic decision tree from calculated stats
    # to update the override parameters in place
    limits = {
        "mask_frac": 0.1,
        "max_offset_low_mask": 0.02,
        "max_offset_high_mask": 0.07,
        "max_slope_low_offset": 1e-5,
        "max_slope_high_offset": 1e-5,
        "slope_ratio": 3.0,
        "channel_sigma": 0.02,
    }
    _image_decision_tree(
        parameters, limits, mask_frac, max_offset, max_slope, slope_ratio, channel_sigma
    )

    return parameters
