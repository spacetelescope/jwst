import logging

import numpy as np
from photutils.aperture import ApertureStats, CircularAperture
from scipy.ndimage import median_filter

log = logging.getLogger(__name__)


__all__ = ["center_from_ta_image", "NoFinitePixelsError"]


class NoFinitePixelsError(Exception):
    """Custom exception raised when no finite pixels are found in the TA image."""

    pass


def center_from_ta_image(ta_image, ref_center, subarray_origin=(1, 1)):
    """
    Determine the center of a point source from a TA image.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.
    ref_center : tuple of float
        (x_ref, y_ref) reference center position in full-frame detector coordinates.
    subarray_origin : tuple of int, optional
        (xstart, ystart) 1-indexed origin of the subarray in full-frame coordinates.
        Default is (1, 1) for full frame.

    Returns
    -------
    x_center, y_center : float
        Fitted x, y center position in full-frame detector coordinates.
    """
    # Transform reference center from full-frame to subarray coordinates
    # FITS convention: xstart, ystart are 1-indexed
    # Python/array convention: 0-indexed
    # ref_center is in 0-indexed detector coordinates
    ref_center_subarray = (
        ref_center[0] - (subarray_origin[0] - 1),
        ref_center[1] - (subarray_origin[1] - 1),
    )

    log.info(
        f"Reference center (0-indexed): full-frame=({ref_center[0]:.2f}, {ref_center[1]:.2f}), "
        f"subarray=({ref_center_subarray[0]:.2f}, {ref_center_subarray[1]:.2f})"
    )

    # Create cutout around reference center (in subarray coordinates)
    cutout, cutout_origin = _cutout_center(ta_image, ref_center_subarray)

    if np.sum(~np.isnan(cutout)) < 10:
        raise NoFinitePixelsError(
            "Not enough finite pixels in the cutout for centroid calculation."
        )
    x_center_cutout, y_center_cutout = _fit_centroid(cutout)

    # Transform back to subarray coordinates
    x_center_subarray = x_center_cutout + cutout_origin[0]
    y_center_subarray = y_center_cutout + cutout_origin[1]

    # Transform from subarray to full-frame detector coordinates
    x_center = x_center_subarray + (subarray_origin[0] - 1)
    y_center = y_center_subarray + (subarray_origin[1] - 1)

    log.info(
        f"Fitted center (0-indexed): subarray=({x_center_subarray:.2f}, {y_center_subarray:.2f}), "
        f"full-frame=({x_center:.2f}, {y_center:.2f})"
    )

    return x_center, y_center


def _fit_centroid(ta_image):
    """
    Compute the centroid of the target acquisition image using aperture photometry.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.

    Returns
    -------
    x_center, y_center : float
        Centroid x, y position.
    """
    # Initial guess: spatial median filter to remove hot pixels if present, then find max
    filtered_image = median_filter(ta_image, size=5)
    y_guess, x_guess = np.unravel_index(np.nanargmax(filtered_image), ta_image.shape)
    phot_aper = CircularAperture([x_guess, y_guess], r=10)

    # in-fill NaNs using filtered image for aperture stats calculation
    # aper_stats.centroid cannot handle NaNs
    ta_image_filled = np.where(np.isfinite(ta_image), ta_image, filtered_image)
    aper_stats = ApertureStats(ta_image_filled, phot_aper)
    x_center, y_center = aper_stats.centroid
    return x_center, y_center


def _cutout_center(image, center, size=40):
    """
    Cut out a small square region from an image centered on a reference position.

    Parameters
    ----------
    image : ndarray
        2D image array.
    center : tuple of float
        (x_center, y_center) position for the center of the cutout.
    size : int, optional
        Size of the square cutout in pixels. Default is 20.

    Returns
    -------
    cutout : ndarray
        Square cutout of the image.
    cutout_origin : tuple of int
        (x_origin, y_origin) position of the lower-left corner of the cutout
        in the original image coordinates.
    """
    x_center, y_center = center
    ny, nx = image.shape

    # Convert center to integer pixel for cutout boundaries
    x_center_int = int(np.round(x_center))
    y_center_int = int(np.round(y_center))

    # Calculate half-size
    half_size = size // 2

    # Calculate cutout boundaries
    x_min = max(0, x_center_int - half_size)
    x_max = min(nx, x_center_int + half_size)
    y_min = max(0, y_center_int - half_size)
    y_max = min(ny, y_center_int + half_size)

    # Extract cutout
    cutout = image[y_min:y_max, x_min:x_max]

    # Store the origin for coordinate transformation
    cutout_origin = (x_min, y_min)

    return cutout, cutout_origin
