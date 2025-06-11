"""Subtract a template background reference file from NIRISS SOSS data."""

import numpy as np
from scipy import ndimage
import logging
import warnings

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)

# These parameters were hard-coded in the algorithm derived by NIRISS - future updates
# could consider exposing them as step-level parameters.
BKG_DISCON_COLUMNS = [685, 725]  # The columns to check for a discontinuity in background level.
BACKGROUND_MASK_CUTOFF = 950  # Drop all pixels with col > cutoff for template matching of bkgd.
SUBSTRIP96_ROWSTART = 10  # Determine where to place SUBSTRIP96 slice wrt. SUBSTRIP256 array.


def find_discontinuity(image):
    """
    Determine the location of a discontinuity in the background levels.

    This function applies a gaussian smoothing filter to the image, then
    searches for a discontinuity by locating the maximum gradient in the
    smoothed image. The search is restricted to columns where the
    discontinuity is expected to occur, set by the BKG_DISCON_COLUMNS
    parameter.

    Parameters
    ----------
    image : float32 ndarray
        The background template with a discontinuity.

    Returns
    -------
    float32 ndarray
        The derived location of the discontinuity in [x, y].
    """
    smoothed = ndimage.gaussian_filter(image, 3)
    gradient = np.gradient(smoothed, axis=1)
    x_discontinuity = (
        np.argmax(
            gradient[:, BKG_DISCON_COLUMNS[0] : BKG_DISCON_COLUMNS[1]],
            axis=1,
        )
        + BKG_DISCON_COLUMNS[0]
    )

    return np.column_stack([x_discontinuity, np.arange(len(image))])


def generate_background_masks(background, n_repeats, for_fitting):
    """
    Generate background masks for the two distinct background regions.

    Using the discontinuity in the background template as a separation,
    build selection masks which flag the left and right background
    regions.

    The mask right of the discontinuity is also truncated rightward of
    a cutoff value set by BACKGROUND_MASK_CUTOFF - the NIRISS team
    found that template matching using pixels on the right side of the
    detector led to poorer fits due to source contamination.

    Parameters
    ----------
    background : float32 ndarray
        The 2-D background template.
    n_repeats : int
        The number of integrations in the science data, used to
        broadcast the mask into a 3-D array.
    for_fitting : bool
        If true, the right_mask will be truncated to only use pixels
        left of the BACKGROUND_MASK_CUTOFF value, by default column 950.
        If false, every pixel will belong to either left_mask or
        right_mask.

    Returns
    -------
    bool ndarray, bool ndarray
        The left and right background masks.
    """
    # First find discontinuity in background template
    discontinuity = find_discontinuity(background)

    # Use discontinuity to split base mask into left and right segments
    __, xx = np.indices(background.shape)

    mask_right = xx >= discontinuity[:, 0].reshape(-1, 1)
    mask_left = ~mask_right

    # Remove background pixels rightward of cutoff from mask for template-fitting purposes
    if for_fitting:
        mask_right[xx > BACKGROUND_MASK_CUTOFF] = False

    # Broadcast mask shape to match data integrations
    mask_left = np.repeat(mask_left[np.newaxis, :, :], n_repeats, axis=0)
    mask_right = np.repeat(mask_right[np.newaxis, :, :], n_repeats, axis=0)

    return mask_left, mask_right


def _rms_error(residuals):
    return np.sqrt(np.square(residuals).mean())


def subtract_soss_bkg(
    input_model,
    bkg_name,
    soss_source_percentile,
    soss_bkg_percentile,
):
    """
    Subtract a scaled template reference background from input SOSS data.

    Derive a separate scaling factor left and right of the background discontinuity;
    select the best template from the library of templates in the reference file
    using the root-mean-square of the residuals.

    Parameters
    ----------
    input_model : CubeModel or ImageModel
        The science data, typically multi-integration CubeModel but
        possibly an ImageModel.
    bkg_name : str
        The name of the background reference file.
    soss_source_percentile : float32
        The threshold percentile used as a cutoff - all pixels above the threshold
        are deemed source and are not used for background template matching.
    soss_bkg_percentile : list of float32
        A 2-member list describing the lower and upper limits of the background
        flux percentiles to use for calculation of the scaling factor.

    Returns
    -------
    CubeModel or ImageModel
        The background-subtracted science datamodel.
    """
    # Load background reference file into datamodel
    bkg_refmodel = datamodels.SossBkgModel(bkg_name)

    # Deal with array shape issues - either slice the arrays for SUBSTRIP96 or skip processing.
    if bkg_refmodel.shape[-2] == 256 and input_model.data.shape[-2] == 96:
        # Slice reference model for SUBSTRIP96 data
        bkg_refmodel.data = bkg_refmodel.data[:, SUBSTRIP96_ROWSTART : SUBSTRIP96_ROWSTART + 96, :]
    elif bkg_refmodel.shape[-2] != input_model.data.shape[-2]:
        log.warning(
            "Background reference file shape does not match SOSS subarray. "
            "Skipping subtraction of reference background."
        )
        return None

    # Generate exclusion mask from data array flux threshold, then OR in the DNU dq plane.
    data_mask = input_model.data >= np.nanpercentile(input_model.data, soss_source_percentile)
    data_mask |= input_model.dq & 1 > 0

    # Most SOSS data will be multi-integration - but if input data array is 2-D, cast into
    # 3-D array for ease of computation
    if len(input_model.data.shape) < 3:
        data = input_model.data[np.newaxis, :, :]
        data_mask = data_mask[np.newaxis, :, :]
    else:
        data = input_model.data

    # Initialize best-fit RMSE to large number
    best_rmse = np.inf
    best_fit_template_idx = -1
    best_scales = np.array([np.nan, np.nan])

    # Iterate over template backgrounds in reference file to find best-fit template
    for t_idx, template in enumerate(bkg_refmodel.data):
        # Generate masks for template matching
        masks = generate_background_masks(template, data.shape[-3], for_fitting=True)

        template = np.repeat(template[np.newaxis, :, :], data.shape[-3], axis=0)

        # Initialize some quantities that we'll check for each template
        scales = np.array([np.nan, np.nan])
        rmse = np.array([np.inf, np.inf])

        # upper and lower range of ratios then find median
        for i, mask in enumerate(masks):
            mask[data_mask] = False
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "divide by zero")
                ratio = data[mask] / template[mask]

            q1, q2 = np.nanpercentile(ratio, [soss_bkg_percentile[0], soss_bkg_percentile[1]])
            valid_pixels = (ratio > q1) & (ratio < q2)
            scales[i] = np.nanmedian(ratio[valid_pixels])
            rmse[i] = _rms_error(data[mask] - (template[mask] * scales[i]))

        if np.sum(rmse) < best_rmse:
            best_rmse = np.sum(rmse)
            best_fit_template_idx = t_idx
            best_scales = scales

    if best_fit_template_idx < 0:
        log.warning("Template matching failed for background subtraction. Skipping bkg_subtract.")
        return None

    # With template determined, regenerate masks for separate scaling
    masks = generate_background_masks(
        bkg_refmodel.data[best_fit_template_idx], data.shape[-3], for_fitting=False
    )

    # Special behavior for ImageModels
    if len(input_model.data.shape) < 3:
        masks = [mask[0] for mask in masks]
        template = bkg_refmodel.data[best_fit_template_idx]
    else:
        template = np.repeat(
            bkg_refmodel.data[best_fit_template_idx][np.newaxis, :, :], data.shape[-3], axis=0
        )
    for i, mask in enumerate(masks):
        input_model.data -= template * mask * best_scales[i]

    return input_model
