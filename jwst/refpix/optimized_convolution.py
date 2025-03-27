#
# Module for using the Simple Improved Reference Subtraction (SIRS) algorithm
# to improve the 1/f noise, only for full frame non-IRS2 NIR data
#

import logging
import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_kernels(sirs_kernel_model, detector, gaussmooth, halfwidth):
    """
    Make convolution kernels from Fourier coefficients in the reference file.

    Parameters
    ----------
    sirs_kernel_model : `~jwst.datamodels.SIRSKernelModel`
        Data model containing the Fourier coefficients from the reference files for
        Simple Improved Reference Subtraction (SIRS)
    detector : str
        Name of the detector of the input data
    gaussmooth : float
        Width of Gaussian smoothing kernel to use as a low-pass filter on
        reference file's coefficients
    halfwidth : int
        Half-width of convolution kernel to build from reference file's coefficients

    Returns
    -------
    kernels: list
        List of kernels appropriate for convolution with the left and right reference pixels.
    """
    gamma, zeta = get_conv_kernel_coeffs(sirs_kernel_model, detector)
    if gamma is None or zeta is None:
        log.info(f"Optimized convolution kernel coefficients NOT found for detector {detector}")
        return None

    kernels_left = []
    kernels_right = []
    for chan in range(gamma.shape[0]):
        n = len(gamma[chan]) - 1
        kernel_left = np.fft.fftshift(np.fft.irfft(gamma[chan]))[n - halfwidth : n + halfwidth + 1]
        kernel_right = np.fft.fftshift(np.fft.irfft(zeta[chan]))[n - halfwidth : n + halfwidth + 1]

        x = np.arange(-halfwidth, halfwidth + 1)
        window = np.exp(-(x**2) / (2 * gaussmooth**2))
        window /= np.sum(window)

        kernel_right = np.convolve(kernel_right, window, mode="same")
        kernel_left = np.convolve(kernel_left, window, mode="same")

        kernels_right += [kernel_right]
        kernels_left += [kernel_left]

    return [kernels_left, kernels_right]


def get_conv_kernel_coeffs(sirs_kernel_model, detector):
    """
    Get the convolution kernels coefficients from the reference file.

    Parameters
    ----------
    sirs_kernel_model : `~jwst.datamodels.SIRSKernelModel`
        Data model containing the Fourier coefficients from the reference files for
        Simple Improved Reference Subtraction (SIRS)
    detector : str
        Name of the detector of the input data

    Returns
    -------
    gamma: numpy array
        Fourier coefficients
    zeta: numpy array
        Fourier coefficients
    """
    mdl_dict = sirs_kernel_model.to_flat_dict()
    gamma, zeta = None, None
    for item in mdl_dict:
        det = item.split(sep=".")[0]
        if detector.lower() == det.lower():
            arr_name = item.split(sep=".")[1]
            if arr_name == "gamma":
                gamma = np.array(mdl_dict[item])
            elif arr_name == "zeta":
                zeta = np.array(mdl_dict[item])
        if gamma is not None and zeta is not None:
            break
    return gamma, zeta


def apply_conv_kernel(data, kernels, sigreject=4.0):
    """
    Apply the convolution kernel.

    Parameters
    ----------
    data : 2-D numpy array
        Data to be corrected
    kernels : list
        List containing the left and right kernels
    sigreject : float
        Number of sigmas to reject as outliers

    Returns
    -------
    data : 2-D numpy array
        Data model with convolution
    """
    data = data.astype(float)
    npix = data.shape[-1]

    kernels_l, kernels_r = kernels
    nchan = len(kernels_l)

    l = data[:, :4]
    r = data[:, -4:]

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
    median_absolute_deviation = 1.48
    medl = np.median(l)
    sigl = median_absolute_deviation * np.median(np.abs(l - medl))
    medr = np.median(r)
    sigr = median_absolute_deviation * np.median(np.abs(r - medr))

    # nl and nr are the number of good reference pixels in the left and right
    # channel in each row. These will be used in lieu of replacing the values
    # of those pixels directly.
    goodl = 1 * (np.abs(l - medl) <= sigreject * sigl)
    nl = np.sum(goodl, axis=1)
    goodr = 1 * (np.abs(r - medr) <= sigreject * sigr)
    nr = np.sum(goodr, axis=1)

    # Average of the left and right channels, replacing masked pixels with zeros.
    # Appropriate normalization factors will be computed later.
    l = np.sum(l * goodl, axis=1) / 4
    r = np.sum(r * goodr, axis=1) / 4
    for chan in range(nchan):
        kernel_l = kernels_l[chan]
        kernel_r = kernels_r[chan]

        # Compute normalizations so that we don't have to directly
        # replace the values of flagged/masked reference pixels.
        norm_l = np.convolve(np.ones(nl.shape), kernel_l, mode="same")
        norm_l /= np.convolve(nl / 4, kernel_l, mode="same")
        norm_r = np.convolve(np.ones(nr.shape), kernel_r, mode="same")
        norm_r /= np.convolve(nr / 4, kernel_r, mode="same")
        template = np.convolve(l, kernel_l, mode="same") * norm_l
        template += np.convolve(r, kernel_r, mode="same") * norm_r
        data[:, chan * npix // 4 : (chan + 1) * npix // 4] -= template[:, np.newaxis]

    log.debug("Optimized convolution kernel applied")
    return data
