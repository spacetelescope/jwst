import numpy as np
import logging
from astropy.stats import SigmaClip

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_profile_mask(ref_2d_profile, threshold=1e-3):
    """Build a mask of the trace based on the 2D profile reference file.

    Parameters
    ----------
    ref_2d_profile : array[float]
        The 2d trace profile reference.
    threshold : float
        Threshold value for excluding pixels based on ref_2d_profile.

    Returns
    -------
    bkg_mask : array[bool]
        Pixel mask in the trace based on the 2d profile reference file.
    """

    bkg_mask = (ref_2d_profile > threshold)

    return bkg_mask


def aperture_mask(xref, yref, halfwidth, shape):
    """Build a mask of the trace based on the trace positions.

    Parameters
    ----------
    xref : array[float]
        The reference x-positions.
    yref : array[float]
        The reference y-positions.
    halfwidth : float
        Size of the aperture mask used when extracting the trace
        positions from the data.
    shape : Tuple(int, int)
        The shape of the array to be masked.

    Returns
    -------
    aper_mask : array[bool]
        Pixel mask in the trace based on the given trace positions.
    """

    # Create a coordinate grid.
    x = np.arange(shape[1])
    y = np.arange(shape[0])
    xx, yy = np.meshgrid(x, y)

    # Interpolate the trace positions onto the grid.
    sort = np.argsort(xref)
    ytrace = np.interp(x, xref[sort], yref[sort])

    # Compute the aperture mask.
    aper_mask = np.abs(yy - ytrace) > halfwidth

    return aper_mask


def soss_background(scidata, scimask, bkg_mask=None):
    """Compute a columnwise background for a SOSS observation.

    Parameters
    ----------
    scidata : array[float]
        The image of the SOSS trace.
    scimask : array[bool]
        Boolean mask of pixels to be excluded.
    bkg_mask : array[bool]
        Boolean mask of pixels to be excluded because they are in
        the trace, typically constructed with make_profile_mask.

    Returns
    -------
    scidata_bkg : array[float]
        Background-subtracted image
    col_bkg : array[float]
        Column-wise background values
    npix_bkg : array[float]
        Number of pixels used to calculate each column value in col_bkg
    """

    # Check the validity of the input.
    data_shape = scidata.shape

    if scimask.shape != data_shape:
        msg = 'scidata and scimask must have the same shape.'
        log.critical(msg)
        raise ValueError(msg)

    if bkg_mask is not None:
        if bkg_mask.shape != data_shape:
            msg = 'scidata and bkg_mask must have the same shape.'
            log.critical(msg)
            raise ValueError(msg)

    # Combine the masks and create a masked array.
    if bkg_mask is not None:
        mask = scimask | bkg_mask
    else:
        mask = scimask

    scidata_masked = np.ma.array(scidata, mask=mask)

    # Mask additional pixels using sigma-clipping.
    sigclip = SigmaClip(sigma=3, maxiters=None, cenfunc='mean')
    scidata_clipped = sigclip(scidata_masked, axis=0)

    # Compute the mean for each column and record the number of pixels used.
    col_bkg = scidata_clipped.mean(axis=0)
    col_bkg = np.where(np.all(scidata_clipped.mask, axis=0), 0., col_bkg)
    npix_bkg = (~scidata_clipped.mask).sum(axis=0)

    # Background subtract the science data.
    scidata_bkg = scidata - col_bkg

    return scidata_bkg, col_bkg, npix_bkg


def make_background_mask(deepstack, width=28):
    """Build a mask of the pixels considered to contain the majority of the
    flux, and should therefore not be used to compute the background.

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
    :rtype: array[bool]
    """

    # Get the dimensions of the input image.
    nrows, ncols = np.shape(deepstack)

    # Set the appropriate quantile for masking based on the subarray size.
    if nrows == 96:  # SUBSTRIP96
        quantile = 100 * (1 - width / 96)  # Mask 1 order worth of pixels.
    elif nrows == 256:  # SUBSTRIP256
        quantile = 100 * (1 - 2 * width / 256)  # Mask 2 orders worth of pixels.
    elif nrows == 2048:  # FULL
        quantile = 100 * (1 - 2 * width / 2048)  # Mask 2 orders worth of pixels.
    else:
        msg = (f'Unexpected image dimensions, expected nrows = 96, 256 or 2048, '
               f'got nrows = {nrows}.')
        log.critical(msg)
        raise ValueError(msg)

    # Find the threshold value associated with the quantile.
    threshold = np.nanpercentile(deepstack, quantile)

    # Mask pixels above the threshold value.
    with np.errstate(invalid='ignore'):
        bkg_mask = (deepstack > threshold) | ~np.isfinite(deepstack)

    return bkg_mask


def soss_oneoverf_correction(scidata, scimask, deepstack, bkg_mask=None,
                             zero_bias=False):
    """Compute a column-wise correction to the 1/f noise on the difference image
    of an individual SOSS integration (i.e. an individual integration - a deep
    image of the same observation).

    Parameters
    ----------
    scidata : array[float]
        Image of the SOSS trace.
    scimask : array[boo]
        Boolean mask of pixels to be excluded based on the DQ values.
    deepstack : array[float]
        Deep image of the trace constructed by combining
        individual integrations of the observation.
    bkg_mask : array[bool]
        Boolean mask of pixels to be excluded because they are in the trace,
        typically constructed with make_profile_mask.
    zero_bias : bool
        If True, the corrections to individual columns will be
        adjusted so that their mean is zero.

    Returns
    -------
    scidata_cor : array[float]
        The 1/f-corrected image
    col_cor : array[float]
        The column-wise correction values
    npix_cor : array[float]
        Number of pixels used in each column
    bias : float
        Net change to the image, if zero_bias was False
    """

    # Check the validity of the input.
    data_shape = scidata.shape

    if scimask.shape != data_shape:
        msg = 'scidata and scimask must have the same shape.'
        log.critical(msg)
        raise ValueError(msg)

    if deepstack.shape != data_shape:
        msg = 'scidata and deepstack must have the same shape.'
        log.critical(msg)
        raise ValueError(msg)

    if bkg_mask is not None:

        if bkg_mask.shape != data_shape:
            msg = 'scidata and bkg_mask must have the same shape.'
            log.critical(msg)
            raise ValueError(msg)

    # Subtract the deep stack from the image.
    diffimage = scidata - deepstack

    # Combine the masks and create a masked array.
    mask = scimask | ~np.isfinite(deepstack)

    if bkg_mask is not None:
        mask = mask | bkg_mask

    diffimage_masked = np.ma.array(diffimage, mask=mask)

    # Mask additional pixels using sigma-clipping.
    sigclip = SigmaClip(sigma=3, maxiters=None, cenfunc='mean')
    diffimage_clipped = sigclip(diffimage_masked, axis=0)

    # Compute the mean for each column and record the number of pixels used.
    col_cor = diffimage_clipped.mean(axis=0)
    npix_cor = (~diffimage_clipped.mask).sum(axis=0)

    # Compute the net change to the image.
    bias = np.nanmean(col_cor)

    # Set the net bias to zero.
    if zero_bias:
        col_cor = col_cor - bias

    # Apply the 1/f correction to the image.
    scidata_cor = scidata - col_cor

    return scidata_cor, col_cor, npix_cor, bias
