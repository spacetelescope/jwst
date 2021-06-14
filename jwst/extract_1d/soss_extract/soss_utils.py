#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


def zero_roll(a, shift):
    """Like np.roll but the wrapped around part is set to zero.
    Only works along the first axis of the array.

    :param a: The input array.
    :param shift: The number of rows to shift by.

    :type a: array[any]
    :type shift: int

    :returns: result - the array with the rows shifted.
    :rtype: array[any]
    """

    result = np.zeros_like(a)
    if shift >= 0:
        result[shift:] = a[:-shift]
    else:
        result[:shift] = a[-shift:]

    return result


def robust_polyfit(x, y, order, maxiter=5, nstd=3.):
    """Perform a robust polynomial fit.

    :param x: x data to fit.
    :param y: y data to fit.
    :param order: polynomial order to use.
    :param maxiter: number of iterations for rejecting outliers.
    :param nstd: number of standard deviations to use when rejecting outliers.

    :type x: array[float]
    :type y: array[float]
    :type order: int
    :type maxiter: int
    :type nstd: float

    :returns: param - best-fit polynomial parameters.
    :rtype: array[float]
    """

    mask = np.ones_like(x, dtype='bool')
    for niter in range(maxiter):

        # Fit the data and evaluate the best-fit model.
        param = np.polyfit(x[mask], y[mask], order)
        yfit = np.polyval(param, x)

        # Compute residuals and mask ouliers.
        res = y - yfit
        stddev = np.std(res)
        mask = np.abs(res) <= nstd*stddev

    return param


def get_image_dim(image, header=None, verbose=False):
    """Determine the properties of the image array.

    :param image: A 2D image of the detector.
    :param header: The header from one of the SOSS reference files.
    :param verbose: If set True some diagnostic plots will be made.

    :type image: array[float]
    :type header: astropy.io.fits.Header
    :type verbose: bool

    :returns:
    dimx, dimy
        The dimensions of the stack array.
    xos, yos
        The oversampling factors of the stack array.
    xnative, ynative
        The dimensions of the stack image, in native pixels.
    padding
        Amount of padding around the image, in native pixels.
    refpix_mask
        Boolean array indicating which pixels are lightsensitive (True) and which are reference pixels (False).

    :rtype: Tuple(int, int, int, int, int, int, int, array[bool])
    """

    # Dimensions of the subarray.
    dimy, dimx = np.shape(image)

    # If no header was passed we have to check all possible sizes.
    if header is None:

        # Initialize padding to zero in this case because it is not a reference file.
        padding = 0

        # Assume the stack is a valid SOSS subarray.
        # FULL: 2048x2048 or 2040x2040 (working pixels) or multiple if oversampled.
        # SUBSTRIP96: 2048x96 or 2040x96 (working pixels) or multiple if oversampled.
        # SUBSTRIP256: 2048x256 or 2040x252 (working pixels) or multiple if oversampled.

        # Check if the size of the x-axis is valid.
        if (dimx % 2048) == 0:
            xnative = 2048
            xos = int(dimx // 2048)

        elif (dimx % 2040) == 0:
            xnative = 2040
            xos = int(dimx // 2040)

        else:
            msg = ('Stack X dimension has unrecognized size of {:}. '
                   'Accepts 2048, 2040 or multiple of.')
            raise ValueError(msg.format(dimx))

        # Check if the y-axis is consistent with the x-axis.
        if np.int(dimy/xos) in [96, 256, 252, 2040, 2048]:
            yos = np.copy(xos)
            ynative = np.int(dimy/yos)

        else:
            msg = ('Stack Y dimension ({:}) is inconsistent with '
                   'stack X dimension ({:}) for acceptable SOSS arrays')
            raise ValueError(msg.format(dimy, dimx))

        # Create a boolean mask indicating which pixels are not reference pixels.
        refpix_mask = np.ones_like(image, dtype='bool')
        if xnative == 2048:
            # Mask out the left and right columns of reference pixels.
            refpix_mask[:, :xos * 4] = False
            refpix_mask[:, -xos * 4:] = False

        if ynative == 2048:
            # Mask out the top and bottom rows of reference pixels.
            refpix_mask[:yos * 4, :] = False
            refpix_mask[-yos * 4:, :] = False

        if ynative == 256:
            # Mask the top rows of reference pixels.
            refpix_mask[-yos * 4:, :] = False

    else:
        # Read the oversampling and padding from the header.
        padding = int(header['PADDING'])
        xos, yos = int(header['OVERSAMP']), int(header['OVERSAMP'])

        # Check that the stack respects its intended format.
        if (dimx/xos - 2*padding) not in [2048]:
            raise ValueError('The header passed is inconsistent with the X dimension of the stack.')
        else:
            xnative = 2048

        if (dimy/yos - 2*padding) not in [96, 256, 2048]:
            raise ValueError('The header passed is inconsistent with the Y dimension of the stack.')
        else:
            ynative = np.int(dimy/yos - 2*padding)

        # The trace file contains no reference pixels so all pixels are good.
        refpix_mask = np.ones_like(image, dtype='bool')

    # If verbose print the output.
    if verbose:
        print('Data dimensions:')
        str_args = dimx, dimy, xos, yos, xnative, ynative
        msg = 'dimx={:}, dimy={:}, xos={:}, yos={:}, xnative={:}, ynative={:}'
        print(msg.format(*str_args))

    return dimx, dimy, xos, yos, xnative, ynative, padding, refpix_mask


def main():
    """Placeholder for potential multiprocessing."""

    return


if __name__ == '__main__':
    main()
