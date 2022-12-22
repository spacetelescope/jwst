import logging
import numpy as np
import warnings

from .soss_utils import robust_polyfit, get_image_dim

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def center_of_mass(column, ypos, halfwidth):
    """Compute a windowed center-of-mass along a column.

    Parameters
    ----------
    column : array[float]
        The column values on which to compute the windowed center of mass.
    ypos : float
        The position along the column to center the window on.
    halfwidth : int
        The half-size of the window in pixels.

    Returns
    --------
    ycom : float
        The center-of-mass of the pixels within the window.
    """

    # Get the column shape and create a corresponding array of positions.
    dimy, = column.shape
    ypix = np.arange(dimy)

    # Find the indices of the window.
    miny = int(np.fmax(np.around(ypos - halfwidth), 0))
    maxy = int(np.fmin(np.around(ypos + halfwidth + 1), dimy))

    # Compute the center of mass on the window.
    with np.errstate(invalid='ignore'):
        ycom = (np.nansum(column[miny:maxy] * ypix[miny:maxy]) /
                np.nansum(column[miny:maxy]))

    return ycom


def get_centroids_com(scidata_bkg, header=None, mask=None, poly_order=11):
    """Determine the x, y coordinates of the trace using a center-of-mass
    analysis. Works for either order if there is no contamination, or for
    order 1 on a detector where the two orders are overlapping.

    Parameters
    ----------
    scidata_bkg : array[float]
        A background subtracted observation.
    header : astropy.io.fits.Header
        The header from one of the SOSS reference files.
    mask : array[bool]
        A boolean array of the same shape as image. Pixels corresponding to
        True values will be masked.
    poly_order : None or int
        Order of the polynomial to fit to the extracted trace positions.

    Returns
    --------
    xtrace : array[float]
        The x coordinates of trace as computed from the best fit polynomial.
    ytrace : array[float]
        The y coordinates of trace as computed from the best fit polynomial.
    param : array[float]
        The best-fit polynomial parameters.
    """

    # If no mask was given use all pixels.
    if mask is None:
        mask = np.zeros_like(scidata_bkg, dtype='bool')

    # Call the script that determines the dimensions of the stack.
    dimx, dimy, xos, yos, xnative, ynative, padding, refpix_mask = get_image_dim(scidata_bkg, header=header)

    # Replace masked pixel values with NaNs.
    scidata_bkg_masked = np.where(mask | ~refpix_mask, np.nan, scidata_bkg)

    # Find centroid - first pass, use all pixels in the column.

    # Normalize each column
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        maxvals = np.nanmax(scidata_bkg_masked, axis=0)
    scidata_norm = scidata_bkg_masked / maxvals

    # Create 2D Array of pixel positions.
    xpix = np.arange(dimx)
    ypix = np.arange(dimy)
    _, ygrid = np.meshgrid(xpix, ypix)

    # CoM analysis to find initial positions using all rows.
    with np.errstate(invalid='ignore'):
        ytrace = np.nansum(scidata_norm * ygrid, axis=0) / np.nansum(scidata_norm, axis=0)
        ytrace = np.where(np.abs(ytrace) == np.inf, np.nan, ytrace)

    # Second pass - use a windowed CoM at the previous position.
    halfwidth = 30 * yos
    for icol in range(dimx):

        ycom = center_of_mass(scidata_norm[:, icol], ytrace[icol], halfwidth)

        # If NaN was returned or centroid is out of bounds, we are done.
        if not np.isfinite(ycom) or (ycom > (ynative - 1) * yos) or (ycom < 0):
            ytrace[icol] = np.nan
            continue

        # If the pixel at the centroid is below the local mean we are likely
        # mid-way between orders and we should shift the window downward to
        # get a reliable centroid for order 1.
        irow = int(np.around(ycom))
        miny = int(np.fmax(np.around(ycom) - halfwidth, 0))
        maxy = int(np.fmin(np.around(ycom) + halfwidth + 1, dimy))
        if scidata_norm[irow, icol] < np.nanmean(scidata_norm[miny:maxy, icol]):
            ycom = center_of_mass(scidata_norm[:, icol], ycom - halfwidth, halfwidth)

        # If the updated position is too close to the array edge, use NaN.
        if not np.isfinite(ycom) or (ycom <= 5 * yos) or (ycom >= (ynative - 6) * yos):
            ytrace[icol] = np.nan
            continue

        # Update the position if the above checks were successful.
        ytrace[icol] = ycom

    # Third pass - fine tuning using a smaller window.
    halfwidth = 16 * yos
    for icol in range(dimx):

        ytrace[icol] = center_of_mass(scidata_norm[:, icol], ytrace[icol],
                                      halfwidth)

    # Fit the y-positions with a polynomial and use result as true y-positions.
    xtrace = np.arange(dimx)
    mask = np.isfinite(ytrace)

    # For padded arrays ignore padding for consistency with real data
    if padding != 0:
        mask = mask & (xtrace >= xos * padding) & (xtrace < (dimx - xos * padding))

    # If no polynomial order was given return the raw measurements.
    if poly_order is None:
        param = []
    else:
        param = robust_polyfit(xtrace[mask], ytrace[mask], poly_order)
        ytrace = np.polyval(param, xtrace)

    return xtrace, ytrace, param
