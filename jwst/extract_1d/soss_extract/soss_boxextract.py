import numpy as np


def get_box_weights(centroid, n_pix, shape, cols=None):
    """ Return the weights of a box aperture given the centroid and the width of
    the box in pixels. All pixels will have the same weights except at the ends
    of the box aperture.

    :param centroid: Position of the centroid (in rows). Same shape as `cols`
    :param n_pix: Width of the extraction box in pixels.
    :param shape: Shape of the output image. (n_row, n_column)
    :param cols: Column indices of good columns Used if the centroid is defined
        for specific columns or a subrange of columns. # TODO not sure the cols argument adds usefull functionality, remove?

    :type centroid: array[float]
    :type n_pix: float
    :type shape: Tuple(int, int)
    :type cols: array[int]

    :returns: weights - An array of pixel weights to use with the box extraction.
    :rtype: array[float]
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

    :param scidata: 2d array of shape (n_row, n_columns)
        scidata
    :param scierr: 2d array of shape (n_row, n_columns)
        uncertainty map
    :param scimask: 2d array, boolean, same shape as data
        masked pixels
    :param box_weights: 2d array, same shape as data
        pre-computed weights for box extraction.
    :param cols: 1d-array, integer
        Which columns to extract

    :type scidata: array[float]
    :type scierr: array[float]
    :type scimask: array[bool]
    :type box_weights: array[float]
    :type cols: array[int]

    :returns: cols, flux, flux_var - The indices of the extracted columns, the
        flux in each column, and the variance of each column.
    :rtype: Tuple(array[int], array[float], array[float])
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

    # Check that all invalid values are masked.
    if not np.isfinite(data[~mask]).all():
        message = 'scidata contains un-masked invalid values.'
        raise ValueError(message)

    if not np.isfinite(error[~mask]).all():
        message = 'scierr contains un-masked invalid values.'
        raise ValueError(message)

    # Set the weights of masked pixels to zero.
    box_weights[mask] = 0.

    # Extract total flux (sum over columns).
    flux = np.nansum(box_weights*data, axis=0)
    npix = np.nansum(box_weights, axis=0)

    # Extract flux error (sum of variances).
    flux_var = np.nansum(box_weights*error**2, axis=0)
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
    err: 2d array[float]
        Uncertainty map of the pixels.
    data: 2d array[float]
        Pixel values.
    pix_to_estim: 2d array[bool]
        Map of the pixels where the uncertainty need to be estimated.
    valid_pix: 2d array[bool]
        Map of valid pixels to be used to find the error empirically.
    Returns
    -------
    err_filled: 2d array[float]
        same as `err`, but the pixel to be estimated are filled with the estimated values.
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


def main():

    return


if __name__ == '__main__':
    main()
