#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from astropy.stats import SigmaClip


def make_profile_mask(ref_2d_profile, threshold=1e-3):
    """Build a mask of the trace based on the 2D profile reference file.

    :param ref_2d_profile: the 2d trace profile reference.
    :param threshold: threshold value for excluding pixels based on
        ref_2d_profile.
    
    :type ref_2d_profile: array[float]
    :type threshold: float

    :returns: bkg_mask - Masks pixels in the trace based on the 2d profile
        reference file.
    :rtype: array[bool]
    """

    bkg_mask = (ref_2d_profile > threshold)

    return bkg_mask


def aperture_mask(xref, yref, halfwidth, shape):
    """Build a mask of the trace based on the trace positions.

    :param xref: The reference x-positions.
    :param yref: The reference y-positions.
    :param halfwidth: Size of the aperture mask used when extracting the trace
        positions from the data.
    :param shape: The shape of the array to be masked.

    :type xref: array[float]:
    :type yref: array[float]:
    :type halfwidth: float
    :type shape: Tuple(int, int)

    :returns: aper_mask - Masks pixels in the trace based on the given trace
        positions.
    :rtype: array[bool]
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

    :param scidata: the image of the SOSS trace.
    :param scimask: a boolean mask of pixls to be excluded.
    :param bkg_mask: a boolean mask of pixels to be excluded because they are in
        the trace, use for example make_profile_mask to construct such a mask.

    :type scidata: array[float]
    :type scimask: array[bool]
    :type bkg_mask: array[bool]

    :returns: scidata_bkg, col_bkg, npix_bkg - The background subtracted image,
        columnwise background values, and number of pixels used in each column.
    :rtype: Tuple(array[float], array[float], array[float])
    """

    # Check the validity of the input.
    data_shape = scidata.shape

    if scimask.shape != data_shape:
        msg = 'scidata and scimask must have the same shape.'
        raise ValueError(msg)

    if bkg_mask is not None:

        if bkg_mask.shape != data_shape:
            msg = 'scidata and bkg_mask must have the same shape.'
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

    :param deepstack: a deep image of the trace constructed by combining
        individual integrations of the observation.
    :param width: the width of the trace is used to set the fraction of pixels
        to exclude with the mask (i.e. width/256 for a SUBSTRIP256 observation).

    :type deepstack: array[float]
    :type width: int

    :returns: bkg_mask - Masks pixels in the trace based on the deepstack or
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
        msg = ('Unexpected image dimensions, expected nrows = 96, 256 or 2048, '
               'got nrows = {}.')
        raise ValueError(msg.format(nrows))

    # Find the threshold value associated with the quantile.
    threshold = np.nanpercentile(deepstack, quantile)

    # Mask pixels above the threshold value.
    with np.errstate(invalid='ignore'):
        bkg_mask = (deepstack > threshold) | ~np.isfinite(deepstack)  # TODO invalid values in deepstack?

    return bkg_mask


def soss_oneoverf_correction(scidata, scimask, deepstack, bkg_mask=None,
                             zero_bias=False):
    """Compute a columnwise correction to the 1/f noise on the difference image
    of an inidividual SOSS integration (i.e. an individual integration - a deep
    image of the same observation).

    :param scidata: the image of the SOSS trace.
    :param scimask: a boolean mask of pixels to be excluded based on the DQ
        values.
    :param deepstack: a deep image of the trace constructed by combining
        individual integrations of the observation.
    :param bkg_mask: a boolean mask of pixels to be excluded because they are in
        the trace, use for example make_background_mask to construct such a mask.
    :param zero_bias: if True the corrections to individual columns will be
        adjusted so that their mean is zero.

    :type scidata: array[float]
    :type scimask: array[bool]
    :type deepstack: array[float]
    :type bkg_mask: array[bool]
    :type zero_bias: bool

    :returns: scidata_cor, col_cor, npix_cor, bias - The 1/f corrected image,
        columnwise correction values, number of pixels used in each column, and
        the net change to the image if zero_bias was False.
    :rtype: Tuple(array[float], array[float], array[float], float)
    """

    # Check the validity of the input.
    data_shape = scidata.shape

    if scimask.shape != data_shape:
        msg = 'scidata and scimask must have the same shape.'
        raise ValueError(msg)

    if deepstack.shape != data_shape:
        msg = 'scidata and deepstack must have the same shape.'
        raise ValueError(msg)

    if bkg_mask is not None:

        if bkg_mask.shape != data_shape:
            msg = 'scidata and bkg_mask must have the same shape.'
            raise ValueError(msg)

    # Subtract the deep stack from the image.
    diffimage = scidata - deepstack

    # Combine the masks and create a masked array.
    mask = scimask | ~np.isfinite(deepstack)  # TODO invalid values in deepstack?

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


def main():
    """Placeholder for potential multiprocessing."""

    return


if __name__ == '__main__':
    main()
