import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def zero_roll(a, shift):
    """Like np.roll but the wrapped around part is set to zero.
    Only works along the first axis of the array.

    Parameters
    ----------
    a : array
        The input array.
    shift : int
        The number of rows to shift by.

    Returns
    -------
    result : array
        The array with the rows shifted.
    """

    result = np.zeros_like(a)
    if shift >= 0:
        result[shift:] = a[:-shift]
    else:
        result[:shift] = a[-shift:]

    return result


def robust_polyfit(x, y, order, maxiter=5, nstd=3.):
    """Perform a robust polynomial fit.

    Parameters
    ----------
    x : array[float]
        x data to fit.
    y : array[float]
        y data to fit.
    order : int
        polynomial order to use.
    maxiter : int, optional
        number of iterations for rejecting outliers.
    nstd : float, optional
        number of standard deviations to use when rejecting outliers.

    Returns
    -------
    param : array[float]
        best-fit polynomial parameters.
    """

    mask = np.ones_like(x, dtype='bool')
    for niter in range(maxiter):

        # Fit the data and evaluate the best-fit model.
        param = np.polyfit(x[mask], y[mask], order)
        yfit = np.polyval(param, x)

        # Compute residuals and mask outliers.
        res = y - yfit
        stddev = np.std(res)
        mask = np.abs(res) <= nstd * stddev

    return param


def get_image_dim(image, header=None):
    """Determine the properties of the image array.

    Parameters
    ----------
    image : array[float]
        A 2D image of the detector.
    header : astropy.io.fits.Header object, optional
        The header from one of the SOSS reference files.

    Returns
    -------
    dimx, dimy : int
        X and Y dimensions of the stack array.
    xos, yos : int
        Oversampling factors in x and y dimensions of the stack array.
    xnative, ynative : int
        Size of stack image x and y dimensions, in native pixels.
    padding : int
        Amount of padding around the image, in native pixels.
    refpix_mask : array[bool]
        Boolean array indicating which pixels are light-sensitive (True)
        and which are reference pixels (False).
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
            log_message = f'Stack X dimension has unrecognized size of {dimx}. Accepts 2048, 2040 or multiple of.'
            log.critical(log_message)
            raise ValueError(log_message)

        # Check if the y-axis is consistent with the x-axis.
        if int(dimy / xos) in [96, 256, 252, 2040, 2048]:
            yos = np.copy(xos)
            ynative = int(dimy / yos)

        else:
            log_message = f'Stack Y dimension ({dimy}) is inconsistent with stack X' \
                          f'dimension ({dimx}) for acceptable SOSS arrays'
            log.critical(log_message)
            raise ValueError(log_message)

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
        if (dimx / xos - 2 * padding) not in [2048]:
            log_message = 'The header passed is inconsistent with the X dimension of the stack.'
            log.critical(log_message)
            raise ValueError(log_message)
        else:
            xnative = 2048

        if (dimy / yos - 2 * padding) not in [96, 256, 2048]:
            log_message = 'The header passed is inconsistent with the Y dimension of the stack.'
            log.critical(log_message)
            raise ValueError(log_message)
        else:
            ynative = int(dimy / yos - 2 * padding)

        # The trace file contains no reference pixels so all pixels are good.
        refpix_mask = np.ones_like(image, dtype='bool')

    log.debug('Data dimensions:')
    log.debug(f'dimx={dimx}, dimy={dimy}, xos={xos}, yos={yos}, xnative={xnative}, ynative={ynative}')

    return dimx, dimy, xos, yos, xnative, ynative, padding, refpix_mask
