"""
A module that provides algorithms for initial estimation of shifts
based on 2D histograms.

"""

import logging
import numpy as np


__all__ = ['estimate_2dhist_shift']


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _xy_2dhist(imgxy, refxy, r):
    # This code replaces the C version (arrxyzero) from carrutils.c
    # It is about 5-8 times slower than the C version.
    dx = np.subtract.outer(imgxy[:, 0], refxy[:, 0]).ravel()
    dy = np.subtract.outer(imgxy[:, 1], refxy[:, 1]).ravel()
    idx = np.where((dx < r + 0.5) & (dx >= -r - 0.5) &
                   (dy < r + 0.5) & (dy >= -r - 0.5))
    h = np.histogram2d(dx[idx], dy[idx], 2 * r + 1,
                       [[-r - 0.5, r + 0.5], [-r - 0.5, r + 0.5]])
    return h[0].T


def estimate_2dhist_shift(imgxy, refxy, searchrad=3.0):
    """ Create a 2D matrix-histogram which contains the delta between each
        XY position and each UV position. Then estimate initial offset
        between catalogs.
    """
    log.info("Computing initial guess for X and Y shifts...")

    # create ZP matrix
    zpmat = _xy_2dhist(imgxy, refxy, r=searchrad)

    nonzeros = np.count_nonzero(zpmat)
    if nonzeros == 0:
        # no matches within search radius. Return (0, 0):
        log.warning("No matches found within a search radius of {:g} pixels."
                    .format(searchrad))
        return 0.0, 0.0

    elif nonzeros == 1:
        # only one non-zero bin:
        yp, xp = np.unravel_index(np.argmax(zpmat), zpmat.shape)
        maxval = zpmat[yp, xp]
        xp -= searchrad
        yp -= searchrad
        log.info("Found initial X and Y shifts of {:.4g}, {:.4g} "
                 "based on a single non-zero bin and {} matches"
                 .format(xp, yp, int(maxval)))
        return xp, yp

    (xp, yp), fit_status, fit_sl = _find_peak(zpmat, peak_fit_box=5,
                                           mask=zpmat > 0)

    if fit_status.startswith('ERROR'):
        log.warning("No valid shift found within a search radius of {:g} "
                    "pixels.".format(searchrad))
        return 0.0, 0.0

    xp -= searchrad
    yp -= searchrad

    if fit_status == 'WARNING:EDGE':
        log.info("Found peak in the 2D histogram lies at the edge of the "
                 "histogram. Try increasing 'searchrad' for improved results.")

    # Attempt to estimate "significance of detection":
    maxval = zpmat.max()
    zpmat_mask = (zpmat > 0) & (zpmat < maxval)

    if np.any(zpmat_mask):
        bkg = zpmat[zpmat_mask].mean()
        sig = maxval / np.sqrt(bkg)

    flux = int(zpmat[fit_sl].sum())
    log.info("Found initial X and Y shifts of {:.4g}, {:.4g} "
             "with significance of {:.4g} and {:d} matches"
             .format(xp, yp, sig, flux))

    return xp, yp


def _find_peak(data, peak_fit_box=5, mask=None):
    """
    Find location of the peak in an array. This is done by fitting a second
    degree 2D polynomial to the data within a `peak_fit_box` and computing the
    location of its maximum. An initial
    estimate of the position of the maximum will be performed by searching
    for the location of the pixel/array element with the maximum value.

    Parameters
    ----------
    data : numpy.ndarray
        2D data.

    peak_fit_box : int, optional
        Size (in pixels) of the box around the initial estimate of the maximum
        to be used for quadratic fitting from which peak location is computed.
        It is assumed that fitting box is a square with sides of length
        given by ``peak_fit_box``.

    mask : numpy.ndarray, optional
        A boolean type `~numpy.ndarray` indicating "good" pixels in image data
        (`True`) and "bad" pixels (`False`). If not provided all pixels
        in `image_data` will be used for fitting.

    Returns
    -------
    coord : tuple of float
        A pair of coordinates of the peak.

    fit_status : str
        Status of the peak search. Currently the following values can be
        returned:

        - ``'SUCCESS'``: Fit was successful and peak is not on the edge of
          the input array;
        - ``'ERROR:NODATA'``: Not enough valid data to perform the fit; The
          returned coordinate is the center of input array;
        - ``'WARNING:EDGE'``: Peak lies on the edge of the input array.
          Returned coordinates are the result of a discreet search;
        - ``'WARNING:BADFIT'``: Performed fid did not find a maximum. Returned
          coordinates are the result of a discreet search;
        - ``'WARNING:CENTER-OF-MASS'``: Returned coordinates are the result
          of a center-of-mass estimate instead of a polynomial fit. This is
          either due to too few points to perform a fit or due to a
          failure of the polynomial fit.

    fit_box : a tuple of two tuples
        A tuple of two tuples of the form ``((x1, x2), (y1, y2))`` that
        indicates pixel ranges used for fitting (these indices can be used
        directly for slicing input data)

    """
    # check arguments:
    data = np.asarray(data, dtype=np.float64)
    ny, nx = data.shape

    # find index of the pixel having maximum value:
    if mask is None:
        jmax, imax = np.unravel_index(np.argmax(data), data.shape)
        coord = (float(imax), float(jmax))

    else:
        j, i = np.indices(data.shape)
        i = i[mask]
        j = j[mask]

        if i.size == 0:
            # no valid data:
            coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
            return coord, 'ERROR:NODATA', np.s_[0:ny-1, 0:nx-1]

        ind = np.argmax(data[mask])
        imax = i[ind]
        jmax = j[ind]
        coord = (float(imax), float(jmax))

    if data[jmax, imax] < 1:
        # no valid data: we need some counts in the histogram bins
        coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
        return coord, 'ERROR:NODATA', np.s_[0:ny-1, 0:nx-1]

    # choose a box around maxval pixel:
    x1 = max(0, imax - peak_fit_box // 2)
    x2 = min(nx, x1 + peak_fit_box)
    y1 = max(0, jmax - peak_fit_box // 2)
    y2 = min(ny, y1 + peak_fit_box)

    # if peak is at the edge of the box, return integer indices of the max:
    if imax == x1 or imax == x2 or jmax == y1 or jmax == y2:
        return (float(imax), float(jmax)), 'WARNING:EDGE', np.s_[y1:y2, x1:x2]

    # expand the box if needed:
    if (x2 - x1) < peak_fit_box:
        if x1 == 0:
            x2 = min(nx, x1 + peak_fit_box)
        if x2 == nx:
            x1 = max(0, x2 - peak_fit_box)

    if (y2 - y1) < peak_fit_box:
        if y1 == 0:
            y2 = min(ny, y1 + peak_fit_box)
        if y2 == ny:
            y1 = max(0, y2 - peak_fit_box)

    if x2 - x1 == 0 or y2 - y1 == 0:
        # not enough data:
        coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
        return coord, 'ERROR:NODATA', np.s_[y1:y2, x1:x2]

    # fit a 2D 2nd degree polynomial to data:
    xi = np.arange(x1, x2)
    yi = np.arange(y1, y2)
    x, y = np.meshgrid(xi, yi)
    x = x.ravel()
    y = y.ravel()
    v = np.vstack((np.ones_like(x), x, y, x*y, x*x, y*y)).T
    d = data[y1:y2, x1:x2].ravel()
    if mask is not None:
        m = mask[y1:y2, x1:x2].ravel()
        v = v[m]
        d = d[m]

    if d.size == 0 or np.max(d) <= 0:
        coord = ((nx - 1.0) / 2.0, (ny - 1.0) / 2.0)
        return coord, 'ERROR:NODATA', np.s_[y1:y2, x1:x2]

    if d.size < 6:
        # we need at least 6 points to fit a 2D quadratic polynomial
        # attempt center-of-mass instead:
        dt = d.sum()
        xc = np.dot(v[:, 1], d) / dt
        yc = np.dot(v[:, 2], d) / dt
        return (xc, yc), 'WARNING:CENTER-OF-MASS', np.s_[y1:y2, x1:x2]

    try:
        c = np.linalg.lstsq(v, d, rcond=None)[0]
    except np.linalg.LinAlgError:
        print("WARNING: Least squares failed!\n{}".format(c))

        # attempt center-of-mass instead:
        dt = d.sum()
        xc = np.dot(v[:, 1], d) / dt
        yc = np.dot(v[:, 2], d) / dt
        return (xc, yc), 'WARNING:CENTER-OF-MASS', np.s_[y1:y2, x1:x2]

    # find maximum of the polynomial:
    _, c10, c01, c11, c20, c02 = c
    det = 4 * c02 * c20 - c11**2
    if det <= 0 or ((c20 > 0.0 and c02 >= 0.0) or (c20 >= 0.0 and c02 > 0.0)):
        # polynomial does not have max. return maximum value in the data:
        return coord, 'WARNING:BADFIT', np.s_[y1:y2, x1:x2]

    xm = (c01 * c11 - 2.0 * c02 * c10) / det
    ym = (c10 * c11 - 2.0 * c01 * c20) / det

    if 0.0 <= xm <= (nx - 1.0) and 0.0 <= ym <= (ny - 1.0):
        coord = (xm, ym)
        fit_status = 'SUCCESS'

    else:
        xm = 0.0 if xm < 0.0 else min(xm, nx - 1.0)
        ym = 0.0 if ym < 0.0 else min(ym, ny - 1.0)
        fit_status = 'WARNING:EDGE'

    return coord, fit_status, np.s_[y1:y2, x1:x2]
