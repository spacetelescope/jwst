#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings

import numpy as np

from .soss_utils import robust_polyfit, get_image_dim

from matplotlib import colors
import matplotlib.pyplot as plt


def _plot_centroid(image, xtrace, ytrace):
    """Overplot the extracted trace positions on the image.

    Parameters
    ----------
    image : array[float]
        A 2D image of the detector.
    xtrace : array[float]
        The x coordinates of the trace to overplot on the image.
    ytrace : array[float]
        The y coordinates of the trace to overplot on the image.
    """

    nrows, ncols = image.shape

    if nrows == ncols:
        aspect = 1
        figsize = ncols/64, nrows/64
    else:
        aspect = 2
        figsize = ncols/64, nrows/32

    plt.figure(figsize=figsize)

    plt.title('Trace Centroids')

    plt.imshow(image, origin='lower', cmap='inferno', norm=colors.LogNorm(), aspect=aspect)
    plt.plot(xtrace, ytrace, lw=2, ls='--', c='black', label='Centroids')

    plt.xlabel('Spectral Pixel', fontsize=14)
    plt.ylabel('Spatial Pixel', fontsize=14)
    plt.legend(fontsize=12)

    plt.xlim(-0.5, ncols - 0.5)
    plt.ylim(-0.5, nrows - 0.5)

    plt.tight_layout()

    plt.show()
    plt.close()

    return


def center_of_mass(column, ypos, halfwidth):
    """Compute a windowed center-of-mass along a column.

    :param column: The column on which to compute the windowed center of mass.
    :param ypos: The position along the column to center the window on.
    :param halfwidth: The half-size of the window in pixels.

    :type column: array[float]
    :type ypos: float
    :type halfwidth: int

    :returns: ycom - the centerof-mass of the pixels withn the window.
    :rtype: float
    """

    # Get the column shape and create a corresponding array of positions.
    dimy, = column.shape
    ypix = np.arange(dimy)

    # Find the indices of the window.
    miny = np.int(np.fmax(np.around(ypos - halfwidth), 0))
    maxy = np.int(np.fmin(np.around(ypos + halfwidth + 1), dimy))

    # Compute the center of mass on the window.
    with np.errstate(invalid='ignore'):
        ycom = np.nansum(column[miny:maxy]*ypix[miny:maxy])/np.nansum(column[miny:maxy])

    return ycom


def get_centroids_com(scidata_bkg, header=None, mask=None, poly_order=11, verbose=False):
    """Determine the x, y coordinates of the trace using a center-of-mass analysis.
    Works for either order if there is no contamination, or for order 1 on a detector
    where the two orders are overlapping.

    :param scidata_bkg: A background subtracted observation.
    :param header: The header from one of the SOSS reference files.
    :param mask: A boolean array of the same shape as image. Pixels corresponding to True values will be masked.
    :param poly_order: Order of the polynomial to fit to the extracted trace positions.
    :param verbose: If set True some diagnostic plots will be made.

    :type scidata_bkg: array[float]
    :type header: astropy.io.fits.Header
    :type mask: array[bool]
    :type poly_order: None or int
    :type verbose: bool

    :returns: xtrace, ytrace, param - The x, y coordinates of trace as computed from the best fit polynomial
    and the best-fit polynomial parameters.
    :rtype: Tuple(array[float], array[float], array[float])
    """

    # If no mask was given use all pixels.
    if mask is None:
        mask = np.zeros_like(scidata_bkg, dtype='bool')

    # Call the script that determines the dimensions of the stack.
    result = get_image_dim(scidata_bkg, header=header, verbose=verbose)
    dimx, dimy, xos, yos, xnative, ynative, padding, refpix_mask = result

    # Replace masked pixel values with NaNs.
    scidata_bkg_masked = np.where(mask | ~refpix_mask, np.nan, scidata_bkg)

    # Find centroid - first pass, use all pixels in the column.

    # Normalize each column
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        maxvals = np.nanmax(scidata_bkg_masked, axis=0)
    scidata_norm = scidata_bkg_masked/maxvals

    # Create 2D Array of pixel positions.
    xpix = np.arange(dimx)
    ypix = np.arange(dimy)
    _, ygrid = np.meshgrid(xpix, ypix)

    # CoM analysis to find initial positions using all rows.
    with np.errstate(invalid='ignore'):
        ytrace = np.nansum(scidata_norm*ygrid, axis=0)/np.nansum(scidata_norm, axis=0)

    # Second pass - use a windowed CoM at the previous position.
    halfwidth = 30 * yos
    for icol in range(dimx):

        ycom = center_of_mass(scidata_norm[:, icol], ytrace[icol], halfwidth)

        # If NaN was returned we are done.
        if not np.isfinite(ycom):
            ytrace[icol] = np.nan
            continue

        # If the pixel at the centroid is below the local mean we are likely mid-way between orders and
        # we should shift the window downward to get a reliable centroid for order 1.
        irow = np.int(np.around(ycom))
        miny = np.int(np.fmax(np.around(ycom) - halfwidth, 0))
        maxy = np.int(np.fmin(np.around(ycom) + halfwidth + 1, dimy))
        if scidata_norm[irow, icol] < np.nanmean(scidata_norm[miny:maxy, icol]):
            ycom = center_of_mass(scidata_norm[:, icol], ycom - halfwidth, halfwidth)

        # If NaN was returned or the position is too close to the array edge, use NaN.
        if not np.isfinite(ycom) or (ycom <= 5 * yos) or (ycom >= (ynative - 6) * yos):
            ytrace[icol] = np.nan
            continue

        # Update the position if the above checks were succesfull.
        ytrace[icol] = ycom

    # Third pass - fine tuning using a smaller window.
    halfwidth = 16 * yos
    for icol in range(dimx):

        ytrace[icol] = center_of_mass(scidata_norm[:, icol], ytrace[icol], halfwidth)

    # Fit the y-positions with a polynomial and use the result as the true y-positions.
    xtrace = np.arange(dimx)
    mask = np.isfinite(ytrace)

    # For padded arrays ignore padding for consistency with real data
    if padding != 0:
        mask = mask & (xtrace >= xos*padding) & (xtrace < (dimx - xos*padding))

    # If no polynomial order was given return the raw measurements.
    if poly_order is None:
        param = []
    else:
        param = robust_polyfit(xtrace[mask], ytrace[mask], poly_order)
        ytrace = np.polyval(param, xtrace)

    # If verbose visualize the result.
    if verbose is True:
        _plot_centroid(scidata_bkg_masked, xtrace, ytrace)

    return xtrace, ytrace, param


def main():
    """Placeholder for potential multiprocessing."""

    return


if __name__ == '__main__':
    main()
