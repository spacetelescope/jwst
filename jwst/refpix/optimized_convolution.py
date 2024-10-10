#
# Module for using Reference Pixels to improve the 1/f noise, to be
# used be only for non-IRS2 NIR data
#

import logging
import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_kernels(conv_kernel_model, detector, gaussmooth, halfwidth):
    """
    Make convolution kernels reference file's Fourier coefficients.

    Parameters:
    -----------

    conv_kernel_model : `~jwst.datamodels.ConvKernelModel`
        Data model containing the Fourier coefficients from the reference files

    detector : str
        Name of the detector of the input data

    gaussmooth : float
        Width of Gaussian smoothing kernel to use as a low-pass filter on reference file's coefficients

    halfwidth : int
        Half-width of convolution kernel to build from reference file's coefficients

    Returns:
    --------

    kernels: list
        Contains the kernels appropriate for convolution with the left  and right reference pixels

    """

    gamma, zeta = get_conv_kernel_coeffs(conv_kernel_model, detector)
    if gamma is None or zeta is None:
        log.info(f'Optimized convolution kernel coefficients NOT found for detector {detector}')
        return None

    kernel_left = []
    kernel_right = []
    for chan in range(gamma.shape[0]):
        n = len(gamma[chan]) - 1
        kernel_left = np.fft.fftshift(np.fft.irfft(gamma[chan]))[n - halfwidth:n + halfwidth + 1]
        kernel_right = np.fft.fftshift(np.fft.irfft(zeta[chan]))[n - halfwidth:n + halfwidth + 1]

        x = np.arange(-halfwidth, halfwidth + 1)
        window = np.exp(-x ** 2 / (2 * gaussmooth ** 2))
        window /= np.sum(window)

        kernel_right = np.convolve(kernel_right, window, mode='same')
        kernel_left = np.convolve(kernel_left, window, mode='same')

        kernel_right += kernel_right
        kernel_left += kernel_left

    return [kernel_left, kernel_right]


def get_conv_kernel_coeffs(conv_kernel_model, detector):
    """
    Get the convolution kernels coefficients from the reference file

    Parameters:
    -----------

    conv_kernel_model : `~jwst.datamodels.ConvKernelModel`
        Data model containing the Fourier coefficients from the reference files

    detector : str
        Name of the detector of the input data

    Returns:
    --------

    gamma: numpy array
        Fourier coefficients

    zeta: numpy array
        Fourier coefficients
    """
    mdl_dict = conv_kernel_model.to_flat_dict()
    gamma, zeta = None, None
    for item in mdl_dict:
        det = item.split(sep='.')[0]
        if detector.lower() == det.lower():
            arr_name = item.split(sep='.')[1]
            if arr_name == 'gamma':
                gamma = np.array(mdl_dict[item])
            elif arr_name == 'zeta':
                zeta = np.array(mdl_dict[item])
        if gamma is not None and zeta is not None:
            break
    return gamma, zeta


def apply_conv_kernel(data, kernels, zeroim, sigreject=4.0):
    """
    Apply the convolution kernel.

    Parameters:
    -----------

    data : 2-D numpy array
        Data to be corrected

    kernels : list
        List containing the left and right kernels

    zeroim : 2-D numpy array
        First group of first integration, to find outliers

    sigreject: float
        Number of sigmas to reject as outliers

    Returns:
    --------

    data : 2-D numpy array
        Data model with convolution
    """
    data = data.astype(float)
    npix = data.shape[-1]

    kernels_l, kernels_r = kernels
    nchan = len(kernels_l)

    # The subtraction below is needed to flag outliers in the reference pixels.
    zeroim = zeroim.astype(float)
    data -= zeroim

    L = data[:, :4]
    R = data[:, -4:]

    # Find the approximate standard deviations of the reference pixels
    # using an outlier-robust median approach. Mask pixels that differ
    # by more than sigreject sigma from this level.
    # NOTE: The Median Absolute Deviation (MAD) is calculated as the
    #   median of the absolute differences between data values and their
    #   median. For normal distribution MAD is equal to 1.48 times the
    #   standard deviation but is a more robust estimate of the dispersion
    #   of data values.The calculation of MAD is straightforward but
    #   time-consuming, especially if MAD estimates are needed for the
    #   local environment around every pixel of a large image. The
    #   calculation is MAD = np.median(np.abs(x-np.median(x))).
    #   Reference: https://www.interstellarmedium.org/numerical_tools/mad/
    MAD = 1.48
    medL = np.median(L)
    sigL = MAD * np.median(np.abs(L - medL))
    medR = np.median(R)
    sigR = MAD * np.median(np.abs(R - medR))

    # nL and nR are the number of good reference pixels in the left and right
    # channel in each row. These will be used in lieu of replacing the values
    # of those pixels directly.
    goodL = 1 * (np.abs(L - medL) <= sigreject * sigL)
    nL = np.sum(goodL, axis=1)
    goodR = 1 * (np.abs(R - medR) <= sigreject * sigR)
    nR = np.sum(goodR, axis=1)

    # Average of the left and right channels, replacing masked pixels with zeros.
    # Appropriate normalization factors will be computed later.
    L = np.sum(L * goodL, axis=1) / 4
    R = np.sum(R * goodR, axis=1) / 4
    for chan in range(nchan):
        kernel_l = kernels_l[chan]
        kernel_r = kernels_r[chan]

        # Compute normalizations so that we don't have to directly
        # replace the values of flagged/masked reference pixels.
        normL = np.convolve(np.ones(nL.shape), kernel_l, mode='same')
        normL /= np.convolve(nL / 4, kernel_l, mode='same')
        normR = np.convolve(np.ones(nR.shape), kernel_r, mode='same')
        normR /= np.convolve(nR / 4, kernel_r, mode='same')
        template = np.convolve(L, kernel_l, mode='same') * normL
        template += np.convolve(R, kernel_r, mode='same') * normR
        data[:, chan * npix * 3 // 4:(chan + 1) * npix * 3 // 4] -= template[:, np.newaxis]

    data += zeroim
    log.debug('Optimized convolution kernel applied')
    return data

