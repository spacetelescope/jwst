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


def main():

    return


if __name__ == '__main__':
    main()
