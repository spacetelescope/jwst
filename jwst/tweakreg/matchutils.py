"""
A module that provides algorithms for initial estimation of shifts
based on 2D histograms.

"""

import logging
import numpy as np
import stsci.imagestats as imagestats

# LOCAL
try:
    from . import chelp
except ImportError:
    chelp = None
from . import __version__
from . import __vdate__

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def center_of_mass(img, labels=None, index=None):
    """
    Calculate the center of mass of the values of an array at labels.

    Parameters
    ----------
    img : ndarray
        Data from which to calculate center-of-mass.

    Returns
    -------
    centerofmass : tuple, or list of tuples
        Coordinates of centers-of-masses.

    Examples
    --------
    >>> from jwst.tweakreg import matchutils
    >>> a = np.array(([0,0,0,0],
                      [0,1,1,0],
                      [0,1,1,0],
                      [0,1,1,0]))
    >>> matchutils.center_of_mass(a)
    (2.0, 1.5)

    """
    normalizer = img.sum(dtype=np.float64)
    if normalizer == 0.0:
        invnorm = np.nan
    else:
        invnorm = 1.0 / normalizer

    grids = np.ogrid[[slice(0, i) for i in img.shape]]

    results = [(img * grids[d].astype(np.float)).sum(dtype=np.float64) *
               invnorm for d in range(img.ndim)]

    if np.isscalar(results[0]):
        return tuple(results)

    return [tuple(v) for v in np.array(results).T]


def build_xy_zeropoint(imgxy, refxy, searchrad=3.0):
    """ Create a matrix which contains the delta between each XY position and
        each UV position.
    """
    log.info("Computing initial guess for X and Y shifts...")

    if chelp is None:
        raise ImportError('cannot import chelp')

    # run C function to create ZP matrix
    zpmat = chelp.arrxyzero(imgxy.astype(np.float32),
                            refxy.astype(np.float32), searchrad)

    xp, yp, flux, zpqual = find_xy_peak(zpmat, center=(searchrad, searchrad))
    if zpqual is None:
        # try with a lower sigma to detect a peak in a sparse set of sources
        xp, yp, flux, zpqual = find_xy_peak(
            zpmat,
            center=(searchrad, searchrad),
            sigma=1.0
        )

    if zpqual is None:
        log.warning("No valid shift found within a search radius of {:g} "
                    "pixels.".format(searchrad))
    else:
        log.info("Found initial X and Y shifts of {:.4g}, {:.4g} "
                 "with significance of {:.4g} and {} matches"
                 .format(xp, yp, zpqual, flux))

    return xp, yp, flux, zpqual


def find_xy_peak(img, center=None, sigma=3.0):
    """ Find the center of the peak of offsets
    """
    # find level of noise in histogram
    istats = imagestats.ImageStats(img, nclip=1,
                                   fields='stddev,mode,mean,max,min')

    if istats.stddev == 0.0:
        istats = imagestats.ImageStats(img, fields='stddev,mode,mean,max,min')

    imgsum = img.sum()

    # clip out all values below mean+3*sigma from histogram
    imgc = img[:, :].copy()
    imgc[imgc < istats.mode + istats.stddev * sigma] = 0.0

    # identify position of peak
    yp0, xp0 = np.where(imgc == imgc.max())

    # Perform bounds checking on slice from img
    ymin = max(0, int(yp0[0]) - 3)
    ymax = min(img.shape[0], int(yp0[0]) + 4)
    xmin = max(0, int(xp0[0]) - 3)
    xmax = min(img.shape[1], int(xp0[0]) + 4)

    # take sum of at most a 7x7 pixel box around peak
    xp_slice = (slice(ymin, ymax), slice(xmin, xmax))
    yp, xp = center_of_mass(img[xp_slice])

    if np.isnan(xp) or np.isnan(yp):
        xp = 0.0
        yp = 0.0
        flux = 0.0
        zpqual = None

    else:
        xp += xp_slice[1].start
        yp += xp_slice[0].start

        # compute S/N criteria for this peak: flux/sqrt(mean of rest of array)
        flux = imgc[xp_slice].sum()
        delta_size = float(img.size - imgc[xp_slice].size)
        if delta_size == 0:
            delta_size = 1
        delta_flux = float(imgsum - flux)
        if flux > imgc[xp_slice].max():
            delta_flux = flux - imgc[xp_slice].max()
        else:
            delta_flux = flux
        zpqual = flux / np.sqrt(delta_flux / delta_size)
        if np.isnan(zpqual) or np.isinf(zpqual):
            zpqual = None

        if center is not None:
            xp -= center[0]
            yp -= center[1]
        flux = imgc[xp_slice].max()

    return xp, yp, flux, zpqual
