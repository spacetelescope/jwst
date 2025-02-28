import numpy as np
import logging
from astropy.stats import SigmaClip

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def soss_background(scidata, scimask, bkg_mask):
    """
    Compute a columnwise background for a SOSS observation.

    Parameters
    ----------
    scidata : array[float]
        The image of the SOSS trace.
    scimask : array[bool]
        Boolean mask of pixels to be excluded.
    bkg_mask : array[bool]
        Boolean mask of pixels to be excluded because they are in
        the trace, typically constructed with make_background_mask.

    Returns
    -------
    scidata_bkg : array[float]
        Background-subtracted image
    col_bkg : array[float]
        Column-wise background values
    """
    # Check the validity of the input.
    data_shape = scidata.shape

    if (scimask.shape != data_shape) or (bkg_mask.shape != data_shape):
        msg = "scidata, scimask, and bkg_mask must all have the same shape."
        log.critical(msg)
        raise ValueError(msg)

    # Combine the masks and create a masked array.
    mask = scimask | bkg_mask
    scidata_masked = np.ma.array(scidata, mask=mask)

    # Mask additional pixels using sigma-clipping.
    sigclip = SigmaClip(sigma=3, maxiters=None, cenfunc="mean")
    scidata_clipped = sigclip(scidata_masked, axis=0)

    # Compute the mean for each column and record the number of pixels used.
    col_bkg = scidata_clipped.mean(axis=0)
    col_bkg = np.where(np.all(scidata_clipped.mask, axis=0), 0.0, col_bkg)

    # Background subtract the science data.
    scidata_bkg = scidata - col_bkg

    return scidata_bkg, col_bkg


def make_background_mask(deepstack, width):
    """
    Build mask of pixels containing most of the flux.

    These pixels should not be used to compute the background.

    Parameters
    ----------
    deepstack : array[float]
        Deep image of the trace constructed by combining
        individual integrations of the observation.
    width : int
        Width, in pixels, of the trace to exclude with the mask
        (i.e. width/256 for a SUBSTRIP256 observation).

    Returns
    -------
    bkg_mask : array[bool]
        Pixel mask in the trace based on the deepstack or
        non-finite in the image.
    """
    # Get the dimensions of the input image.
    nrows, _ = np.shape(deepstack)

    # Set the appropriate quantile for masking based on the subarray size.
    if nrows == 96:  # SUBSTRIP96
        quantile = 100 * (1 - width / nrows)  # Mask 1 order worth of pixels.
    elif nrows in [256, 2048]:  # SUBSTRIP256, FULL
        quantile = 100 * (1 - 2 * width / nrows)  # Mask 2 orders worth of pixels.
    else:
        msg = f"Unexpected image dimensions, expected nrows = 96, 256 or 2048, got nrows = {nrows}."
        log.critical(msg)
        raise ValueError(msg)

    # Find the threshold value associated with the quantile.
    threshold = np.nanpercentile(deepstack, quantile)

    # Mask pixels above the threshold value.
    with np.errstate(invalid="ignore"):
        return (deepstack > threshold) | ~np.isfinite(deepstack)
