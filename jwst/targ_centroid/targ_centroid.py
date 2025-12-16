import logging

import numpy as np
from photutils.centroids import centroid_2dg

log = logging.getLogger(__name__)


__all__ = ["center_from_ta_image", "NoFinitePixelsError", "BadFitError"]


class NoFinitePixelsError(Exception):
    """Custom exception raised when no finite pixels are found in the TA image."""

    pass


class BadFitError(Exception):
    """Custom exception raised when the model fit does not meet quality criteria."""

    pass


def center_from_ta_image(ta_image, ref_center, subarray_origin=(1, 1)):
    """
    Determine the center of a point source from a TA image.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.
    ref_center : tuple of float
        (x_ref, y_ref) reference center position in subarray coordinates, zero-indexed.
        Dither position comes out of the WCS transform already in subarray coordinates.
    subarray_origin : tuple of int, optional
        (xstart, ystart) 1-indexed origin of the subarray in full-frame coordinates.
        Default is (1, 1) for full frame.

    Returns
    -------
    x_center, y_center : float
        Fitted x, y center position in full-frame detector coordinates.
    """
    log.info("Computing centroid of source in TA verification image.")

    # Create small cutout around reference center (in subarray coordinates)
    cutout, cutout_origin = _cutout_center(ta_image, ref_center, size=20)

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

    log.debug(
        f"Fitted center (0-indexed): subarray=({x_center_subarray:.2f}, {y_center_subarray:.2f}), "
        f"full-frame=({x_center:.2f}, {y_center:.2f})"
    )

    return x_center, y_center


def _fit_centroid(cutout):
    """
    Compute the centroid of the target acquisition image.

    Parameters
    ----------
    cutout : ndarray
        2D target acquisition image data, cut out around the reference position.

    Returns
    -------
    x_center, y_center : float
        Centroid x, y position.
    """
    mask = ~np.isfinite(cutout)
    # Use a 2-D Gaussian fit to find the centroid
    try:
        x_center, y_center = centroid_2dg(cutout, mask=mask)
    except ValueError as e:
        raise BadFitError(
            "2D Gaussian centroid fit failed. Check input data and mask. "
            f"Error from fitter was {type(e).__name__}: {e}"
        ) from None

    return x_center, y_center


def _cutout_center(image, center, size=16):
    """
    Cut out a small square region from an image centered on a reference position.

    Parameters
    ----------
    image : ndarray
        2D image array.
    center : tuple of float
        (x_center, y_center) position for the center of the cutout.
    size : int, optional
        Size of the square cutout in pixels.

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
