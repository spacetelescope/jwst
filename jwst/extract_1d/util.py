import logging

import numpy as np

from gwcs.wcstools import grid_from_bounding_box

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def pixel_area(wcs, shape, is_wfss, verbose=True):
    """Compute the solid angle of a pixel.

    Parameters
    ----------
    wcs : function
        This takes a pair of pixel coordinates (`x`, `y`) and returns
        the right ascension, declination, and wavelength at that pixel
        (or array of pixels).

    shape : tuple
        The shape of the data array.  The bounding box from the wcs will be
        used instead, if it exists.

    is_wfss : bool
        True if the input file contains WFSS data.

    verbose : bool
        If True, log messages.

    Returns
    -------
    pixel_solid_angle : float or None
        The average solid angle of a pixel, in steradians.  None will be
        returned if this value could not be determined.
    """

    if is_wfss:
        if verbose:
            log.warning("The solid angle of a pixel cannot currently "
                        "be determined for WFSS data.")
        return None

    if shape is not None and len(shape) != 2 and len(shape) != 3:
        if verbose:
            log.warning("shape = {}, should be 2-D or 3-D".format(shape))
        return None

    three_d = (len(shape) == 3)

    if three_d:
        grid = np.indices(shape[-2:])
        z = np.zeros(shape[-2:], dtype=np.float64) + shape[0] // 2
        # The arguments are the X, Y, and Z pixel coordinates, and the
        # output arrays will be 2-D.
        try:
            stuff = wcs(grid[1], grid[0], z)
        except TypeError:
            stuff = temp_wcs(wcs, grid[1], grid[0], z)
        (ra, dec) = stuff[0:2]
    else:
        if hasattr(wcs, 'bounding_box') and wcs.bounding_box is not None:
            grid2 = grid_from_bounding_box(wcs.bounding_box)
            # Note that the elements of grid2 are reversed wrt grid.
            try:
                stuff = wcs(grid2[0], grid2[1])
            except TypeError:
                stuff = temp_wcs(wcs, grid2[0], grid2[1])
            (ra, dec) = stuff[0:2]
        else:
            grid = np.indices(shape)
            try:
                stuff = wcs(grid[1], grid[0])
                (ra, dec) = stuff[0:2]
            except TypeError:
                stuff = temp_wcs(wcs, grid[1], grid[0])
            (ra, dec) = stuff[0:2]

    cos_dec = np.cos(dec[0:-1, 0:-1] * np.pi / 180.)
    dra1 = (ra[0:-1, 1:] - ra[0:-1, 0:-1]) * cos_dec
    dra2 = (ra[1:, 0:-1] - ra[0:-1, 0:-1]) * cos_dec
    ddec1 = dec[0:-1, 1:] - dec[0:-1, 0:-1]
    ddec2 = dec[1:, 0:-1] - dec[0:-1, 0:-1]

    cdelt1 = np.sqrt(dra1**2 + ddec1**2)
    cdelt2 = np.sqrt(dra2**2 + ddec2**2)
    cdelt1 = np.nanmean(cdelt1)
    cdelt2 = np.nanmean(cdelt2)
    if verbose:
        log.debug("Pixel dimensions in arcsec:  {:.5g} x {:.5g}"
                  .format(cdelt1 * 3600., cdelt2 * 3600.))

    if cdelt1 == 0. and cdelt2 == 0:
        if verbose:
            log.warning("Pixel solid angle could not be determined")
        pixel_solid_angle = None
    else:
        if cdelt1 < 0.01 * cdelt2:
            if verbose:
                log.warning("Using cdelt2 = {:.4g} arcsec for both dimensions "
                            "of a pixel".format(cdelt2 * 3600.))
            cdelt1 = cdelt2
        elif cdelt2 < 0.01 * cdelt1:
            if verbose:
                log.warning("Using cdelt1 = {:.4g} arcsec for both dimensions "
                            "of a pixel".format(cdelt1 * 3600.))
            cdelt2 = cdelt1
        pixel_solid_angle = (cdelt1 * cdelt2) * (np.pi / 180.)**2

    return pixel_solid_angle


def temp_wcs(wcs, x, y, z=None):
    """Call the wcs function in a loop over pixels.

    Parameters
    ----------
    wcs : function
        This takes a pair of pixel coordinates and returns the right
        ascension, declination, and wavelength at that pixel.

    x, y : ndarray, 2-D
        The arrays of pixel coordinates.

    z : ndarray, or None
        If the data are 3-D, e.g. IFU or CubeModel, this should be an array
        of values for one plane in axis 0 (wavelength for IFU, multiple
        exposures for CubeModel).

    Returns
    -------
    ra, dec, wl : ndarray, 2-D
        The first three elements of the output of the `wcs` function for
        spectral data.
    """

    ra = np.zeros_like(x)
    dec = np.zeros_like(x)
    wl = np.zeros_like(x)
    shape = x.shape
    if z is None:
        for j in range(shape[-2]):
            for i in range(shape[-1]):
                stuff = wcs(x[j, i], y[j, i])
                ra[j, i], dec[j, i], wl[j, i] = stuff[0:3]
    else:
        for j in range(shape[-2]):
            for i in range(shape[-1]):
                stuff = wcs(x[j, i], y[j, i], z[j, i])
                ra[j, i], dec[j, i], wl[j, i] = stuff[0:3]

    return (ra, dec, wl)
