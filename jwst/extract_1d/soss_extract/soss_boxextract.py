import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def get_box_weights(centroid, n_pix, shape, cols=None):
    """ Return the weights of a box aperture given the centroid and the width of
    the box in pixels. All pixels will have the same weights except at the ends
    of the box aperture.

    Parameters
    ----------
    centroid : array[float]
        Position of the centroid (in rows). Same shape as `cols`
    n_pix : float
        Width of the extraction box in pixels.
    shape : Tuple(int, int)
        Shape of the output image. (n_row, n_column)
    cols : array[int]
        Column indices of good columns. Used if the centroid is defined
        for specific columns or a sub-range of columns.

    Returns
    -------
    weights : array[float]
        An array of pixel weights to use with the box extraction.
    """

    nrows, ncols = shape

    # Use all columns if not specified
    if cols is None:
        cols = np.arange(ncols)

    # Row centers of all pixels.
    rows = np.indices((nrows, len(cols)))[0]

    # Pixels that are entierly inside the box are set to one.
    cond = (rows <= (centroid - 0.5 + n_pix / 2))
    cond &= ((centroid + 0.5 - n_pix / 2) <= rows)
    weights = cond.astype(float)

    # Fractional weights at the upper bound.
    cond = (centroid - 0.5 + n_pix / 2) < rows
    cond &= (rows < (centroid + 0.5 + n_pix / 2))
    weights[cond] = (centroid + n_pix / 2 - (rows - 0.5))[cond]

    # Fractional weights at the lower bound.
    cond = (rows < (centroid + 0.5 - n_pix / 2))
    cond &= ((centroid - 0.5 - n_pix / 2) < rows)
    weights[cond] = (rows + 0.5 - (centroid - n_pix / 2))[cond]

    # Return with the specified shape with zeros where the box is not defined.
    out = np.zeros(shape, dtype=float)
    out[:, cols] = weights

    return out


def box_extract(scidata, scierr, scimask, box_weights, cols=None):
    """ Perform a box extraction.

    Parameters
    ----------
    scidata : array[float]
        2d array of science data with shape (n_row, n_columns)
    scierr : array[float]
        2d array of uncertainty map with same shape as scidata
    scimask : array[bool]
        2d boolean array of masked pixels with same shape as scidata
    box_weights : array[float]
        2d array of pre-computed weights for box extraction,
        with same shape as scidata
    cols : array[int]
        1d integer array of column numbers to extract

    Returns
    -------
    cols : array[int]
        Indices of extracted columns
    flux : array[float]
        The flux in each column
    flux_var : array[float]
        The variance of the flux in each column
    """

    nrows, ncols = scidata.shape

    # Use all columns if not specified
    if cols is None:
        cols = np.arange(ncols)

    # Keep only required columns and make a copy.
    data = scidata[:, cols].copy()
    error = scierr[:, cols].copy()
    mask = scimask[:, cols].copy()
    box_weights = box_weights[:, cols].copy()

    # Check that all invalid values in the extraction region are masked.
    extract_region = (box_weights > 0)
    condition = extract_region & ~mask

    if not np.isfinite(data[condition]).all():
        message = 'scidata contains un-masked invalid values.'
        log.critical(message)
        raise ValueError(message)

    if not np.isfinite(error[condition]).all():
        message = 'scierr contains un-masked invalid values.'
        log.critical(message)
        raise ValueError(message)

    # Set all pixels values outside of extraction region to Nan
    # so it will be correctly handle by np.nansum.
    data = np.where(extract_region, data, np.nan)
    error = np.where(extract_region, error, np.nan)

    # Set the weights of masked pixels to zero.
    box_weights[mask] = 0.

    # Extract total flux (sum over columns).
    flux = np.nansum(box_weights * data, axis=0)
    npix = np.nansum(box_weights, axis=0)

    # Extract flux error (sum of variances).
    flux_var = np.nansum(box_weights * error**2, axis=0)
    flux_err = np.sqrt(flux_var)

    # Set empty columns to NaN.
    flux = np.where(npix > 0, flux, np.nan)
    flux_err = np.where(npix > 0, flux_err, np.nan)

    return cols, flux, flux_err, npix


def estim_error_nearest_data(err, data, pix_to_estim, valid_pix):
    """
    Function to estimate pixel error empirically using the corresponding error
    of the nearest pixel value (`data`). Intended to be used in a box extraction
    when the bad pixels are modeled.

    Parameters
    ----------
    err : 2d array[float]
        Uncertainty map of the pixels.
    data : 2d array[float]
        Pixel values.
    pix_to_estim : 2d array[bool]
        Map of the pixels where the uncertainty needs to be estimated.
    valid_pix : 2d array[bool]
        Map of valid pixels to be used to find the error empirically.
    Returns
    -------
    err_filled : 2d array[float]
        same as `err`, but the pixels to be estimated are filled with the estimated values.
    """
    # Tranform to 1d arrays
    data_to_estim = data[pix_to_estim]
    err_valid = err[valid_pix]
    data_valid = data[valid_pix]

    #
    # Use np.searchsorted for efficiency
    #
    # Need to sort the arrays used to find similar values
    idx_sort = np.argsort(data_valid)
    err_valid = err_valid[idx_sort]
    data_valid = data_valid[idx_sort]

    # Searchsorted: gives the position of the nearest higher value,
    # not necessarily the closest value
    idx_higher = np.searchsorted(data_valid, data_to_estim)
    idx_higher = np.clip(idx_higher, 0, err_valid.size - 1)
    # The nearest lower value is given by the preceding index
    idx_lower = np.clip(idx_higher - 1, 0, err_valid.size - 1)

    # Find the best between index around the value (lower and higher index) ...
    idx_around = np.vstack([idx_lower, idx_higher])
    # ... using the one with the smallest error
    distance = np.abs(data_valid[idx_around] - data_to_estim[None, :])
    idx_best_of_2 = np.argmin(distance, axis=0)
    idx_closest = idx_around[idx_best_of_2, np.arange(idx_best_of_2.size)]

    # Get the corresponding error (that's what we want to find!)
    err_estimate = err_valid[idx_closest]

    # Replace estimated values in the ouput error 2d image
    err_out = err.copy()
    err_out[pix_to_estim] = err_estimate

    return err_out
