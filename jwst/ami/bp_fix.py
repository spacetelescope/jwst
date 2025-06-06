"""Pipeline implementation of Jens Kammerer's bp_fix code based on Ireland 2013 algorithm."""

import numpy as np
import logging
import warnings

from copy import deepcopy
from .matrix_dft import matrix_dft
from scipy.ndimage import median_filter

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
logging.captureWarnings(True)

micron = 1.0e-6
filts = ["F277W", "F380M", "F430M", "F480M", "F356W", "F444W"]
# TO DO: get these from some context-controlled place
filtwl_d = {  # pivot wavelengths
    "F277W": 2.776e-6,  # less than Nyquist
    "F380M": 3.828e-6,
    "F430M": 4.286e-6,
    "F480M": 4.817e-6,
    "F356W": 3.595e-6,  # semi-forbidden
    "F444W": 4.435e-6,  # semi-forbidden
}
filthp_d = {  # half power limits
    "F277W": (2.413e-6, 3.142e-6),
    "F380M": (3.726e-6, 3.931e-6),
    "F430M": (4.182e-6, 4.395e-6),
    "F480M": (4.669e-6, 4.971e-6),
    "F356W": (3.141e-6, 4.068e-6),
    "F444W": (3.880e-6, 5.023e-6),
}
WL_OVERSIZEFACTOR = 0.1  # increase filter wl support by this amount to 'oversize' in wl space

DIAM = 6.559348  # / Flat-to-flat distance across pupil in V3 axis
PUPLDIAM = 6.603464  # / Full pupil file size, incl padding.
PUPL_CRC = 6.603464  # / Circumscribing diameter for JWST primary

DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
JUMP_DET = dqflags.pixel["JUMP_DET"]


def create_wavelengths(filtername):
    """
    Extend filter support slightly past half power points.

    Filter transmissions are quasi-rectangular.

    Parameters
    ----------
    filtername : str
        AMI filter name

    Returns
    -------
    tuple
        Tuple of wavelengths. Center, low, and high.
    """
    wl_ctr = filtwl_d[filtername]
    wl_hps = filthp_d[filtername]
    # both positive quantities below - left is lower wl, rite is higher wl
    dleft = (wl_ctr - wl_hps[0]) * (1 + WL_OVERSIZEFACTOR)
    drite = (-wl_ctr + wl_hps[1]) * (1 + WL_OVERSIZEFACTOR)

    return (wl_ctr, wl_ctr - dleft, wl_ctr + drite)


def calc_pupil_support(filtername, sqfov_npix, pxsc_rad, pupil_mask):
    """
    Calculate psf at low, center, and high wavelengths of filter.

    Coadd psfs and perform fft-style transform of image w/ dft.

    Parameters
    ----------
    filtername : str
        AMI filter name
    sqfov_npix : float
        Square field of view in number of pixels
    pxsc_rad : float
        Detector pixel scale in rad/px
    pupil_mask : array
        Pupil mask model (NRM)

    Returns
    -------
    np.array
        Absolute value of FT(image) in filter - the CV Vsq array
    """
    wls = create_wavelengths(filtername)
    log.info(f"      {filtername}: {wls[0] / micron:.3f} to {wls[2] / micron:.3f} micron")
    detimage = np.zeros((sqfov_npix, sqfov_npix), float)
    for wl in wls:
        psf = calcpsf(wl, sqfov_npix, pxsc_rad, pupil_mask)
        detimage += psf

    return transform_image(detimage)


def transform_image(image):
    """
    Take Fourier transform of image.

    Parameters
    ----------
    image : numpy array
        Science image

    Returns
    -------
    np.array
        Absolute value of FT(image)
    """
    ftimage = matrix_dft(
        image,
        image.shape[0],
        image.shape[0],
        centering="ADJUSTABLE",
    )

    return np.abs(ftimage)


def calcpsf(wl, fovnpix, pxsc_rad, pupil_mask):
    """
    Calculate the PSF.

    Parameters
    ----------
    wl : float
        Wavelength (meters)
    fovnpix : float
        Square field of view in number of pixels
    pxsc_rad : float
        Detector pixel scale in rad/px
    pupil_mask : array
        Pupil mask model (NRM)

    Returns
    -------
    image_intensity : numpy array
        Monochromatic unnormalized psf
    """
    reselt = wl / PUPLDIAM  # radian
    nlam_d = fovnpix * pxsc_rad / reselt  # Soummer nlamD FOV in reselts

    image_field = matrix_dft(
        pupil_mask,
        nlam_d,
        fovnpix,
        centering="ADJUSTABLE",
    )

    image_intensity = (image_field * image_field.conj()).real

    return image_intensity


def bad_pixels(data, median_size, median_tres):
    """
    Identify bad pixels by subtracting median-filtered data and searching for outliers.

    Parameters
    ----------
    data : numpy array
        Science data
    median_size : float
        Median filter size (pixels)
    median_tres : float
        Empirically determined threshold

    Returns
    -------
    pxdq : np.ndarray[int]
        Bad pixel mask identified by median filtering
    """
    mfil_data = median_filter(data, size=median_size)
    diff_data = np.abs(data - mfil_data)
    pxdq = diff_data > median_tres * np.median(diff_data)
    pxdq = pxdq.astype("bool")

    log.info(
        f"         Identified {np.sum(pxdq):.0f} bad pixels "
        f"({100.0 * np.sum(pxdq) / np.prod(pxdq.shape):.2f}%)"
    )
    log.info(f"         {np.max(diff_data / np.median(diff_data)):.3f}")

    return pxdq


def fourier_corr(data, pxdq, fmas):
    """
    Compute and apply the bad pixel corrections based on Section 2.5 of Ireland 2013.

    Parameters
    ----------
    data : numpy array
        Science data
    pxdq : numpy array
        Bad pixel mask
    fmas : numpy array
        FT of science data

    Returns
    -------
    data_out : numpy array
        Corrected science data

    References
    ----------
    M. J. Ireland, Phase errors in diffraction-limited imaging: contrast limits
    for sparse aperture masking, Monthly Notices of the Royal Astronomical
    Society, Volume 433, Issue 2, 01 August 2013, Pages 1718–1728,
    https://doi.org/10.1093/mnras/stt859
    """
    # Get the dimensions.
    ww = np.where(pxdq > 0.5)
    ww_ft = np.where(fmas)

    # Compute the B_Z matrix from Section 2.5 of Ireland 2013. This matrix
    # maps the bad pixels onto their Fourier power in the domain Z, which is
    # the complement of the pupil support.
    B_Z = np.zeros((len(ww[0]), len(ww_ft[0]) * 2))  # noqa: N806
    xh = data.shape[0] // 2
    yh = data.shape[1] // 2
    xx, yy = np.meshgrid(
        2.0 * np.pi * np.arange(yh + 1) / data.shape[1],
        2.0 * np.pi * (((np.arange(data.shape[0]) + xh) % data.shape[0]) - xh) / data.shape[0],
    )
    for i in range(len(ww[0])):
        cdft = np.exp(-1j * (ww[0][i] * yy + ww[1][i] * xx))
        B_Z[i, :] = np.append(cdft[ww_ft].real, cdft[ww_ft].imag)

    # Compute the corrections for the bad pixels using the Moore-Penrose pseudo
    # inverse of B_Z (Equation 19 of Ireland 2013).
    B_Z_ct = np.transpose(np.conj(B_Z))  # noqa: N806
    B_Z_mppinv = np.dot(B_Z_ct, np.linalg.inv(np.dot(B_Z, B_Z_ct)))  # noqa: N806

    # Apply the corrections for the bad pixels.
    data_out = deepcopy(data)
    data_out[ww] = 0.0
    data_ft = np.fft.rfft2(data_out)[ww_ft]
    corr = -np.real(np.dot(np.append(data_ft.real, data_ft.imag), B_Z_mppinv))
    data_out[ww] += corr

    return data_out


def fix_bad_pixels(data, pxdq0, filt, pxsc, nrm_model):
    """
    Apply the Fourier bad pixel correction to pixels flagged DO_NOT_USE or JUMP_DET.

    Original code implementation by Jens Kammerer.

    Parameters
    ----------
    data : array
        Cropped science data
    pxdq0 : array
        Cropped DQ array
    filt : str
        AMI filter name
    pxsc : float
        Pixel scale, mas/pixel
    nrm_model : datamodel object
        NRM pupil datamodel

    Returns
    -------
    data : numpy array
        Corrected data
    pxdq : np.ndarray[int]
        Mask of bad pixels, updated if new ones were found
    """
    dq_dnu = pxdq0 & DO_NOT_USE == DO_NOT_USE
    dq_jump = pxdq0 & JUMP_DET == JUMP_DET
    dqmask = dq_dnu | dq_jump

    pxdq = np.where(dqmask, pxdq0, 0)
    nflagged_dnu = np.count_nonzero(pxdq)
    log.info(f"{nflagged_dnu:d} pixels flagged DO_NOT_USE in cropped data")

    # DNU, some other pixels are now NaNs in cal level products.
    # Replace them with 0, then
    # add DO_NOT_USE flags to positions in DQ array so they will be corrected.
    nanidxlist = np.argwhere(np.isnan(data))
    if len(nanidxlist) > 1:
        log.info(f"Identified {len(nanidxlist):d} NaN pixels to correct")
        for idx in nanidxlist:
            data[idx[0], idx[1], idx[2]] = 0
            pxdq0[idx[0], idx[1], idx[2]] += 1  # add DNU flag to each nan pixel

    # These values are taken from the JDox and the SVO Filter Profile
    # Service.
    diam = PUPLDIAM  # m
    gain = 1.61  # e-/ADU
    rdns = 18.32  # e-
    pxsc_rad = (pxsc / 1000) * np.pi / (60 * 60 * 180)

    # These values were determined empirically for NIRISS/AMI and need to be
    # tweaked for any other instrument.
    median_size = 3  # pix
    median_tres = 50.0

    pupil_mask = nrm_model.nrm
    imsz = data.shape
    sh = imsz[-1] // 2  # half size, even
    # Compute field-of-view and Fourier sampling.
    fov = 2 * sh * pxsc / 1000.0  # arcsec
    fsam = filtwl_d[filt] / (fov / 3600.0 / 180.0 * np.pi)  # m/pix
    log.info(f"      FOV = {fov:.1f} arcsec, Fourier sampling = {fsam:.3f} m/pix")

    #
    cvis = calc_pupil_support(filt, 2 * sh, pxsc_rad, pupil_mask)
    cvis /= np.max(cvis)
    fmas = cvis < 1e-3  # 1e-3 seems to be a reasonable threshold
    fmas = np.fft.fftshift(fmas)[:, : 2 * sh // 2 + 1]

    # Compute the pupil mask. This mask defines the region where we are
    # measuring the noise. It looks like 15 lambda/D distance from the PSF
    # is reasonable.
    ramp = np.arange(2 * sh) - 2 * sh // 2
    xx, yy = np.meshgrid(ramp, ramp)
    dist = np.sqrt(xx**2 + yy**2)
    pmas = dist > 9.0 * filtwl_d[filt] / diam * 180.0 / np.pi * 1000.0 * 3600.0 / pxsc

    # Go through all frames.
    for j in range(imsz[0]):
        log.info(f"         Frame {j + 1:.0f} of {imsz[0]:.0f}")

        # Handle odd/even size issues by cropping out the -1th pixel in odd data
        xshape = data.shape[1]
        yshape = data.shape[2]
        if xshape % 2 == 0:
            idx_x = xshape
        elif data.shape[1] % 2 == 1:
            idx_x = xshape - 1
        if yshape % 2 == 0:
            idx_y = yshape
        elif yshape % 2 == 1:
            idx_y = yshape - 1
        data_cut = deepcopy(data[j, :idx_x, :idx_y])
        data_orig = deepcopy(data_cut)
        pxdq_cut = deepcopy(pxdq[j, :idx_x, :idx_y])
        pxdq_cut = pxdq_cut > 0.5
        # Correct the bad pixels. This is an iterative process. After each
        # iteration, we check whether new (residual) bad pixels are
        # identified. If so, we re-compute the corrections. If not, we
        # terminate the iteration.
        for k in range(10):
            # Correct the bad pixels.
            data_cut = fourier_corr(data_cut, pxdq_cut, fmas)
            if k == 0:
                temp = deepcopy(data_cut)

            # Identify residual bad pixels by looking at the high spatial
            # frequency part of the image.
            fmas_data = np.real(np.fft.irfft2(np.fft.rfft2(data_cut) * fmas))

            # Analytically determine the noise (Poisson noise + read noise)
            # and normalize the high spatial frequency part of the image
            # by it, then identify residual bad pixels.
            mfil_data = median_filter(data_cut, size=median_size)
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", category=RuntimeWarning, message="invalid value encountered"
                )
                nois = np.sqrt(mfil_data / gain + rdns**2)
            fmas_data /= nois
            temp = bad_pixels(fmas_data, median_size=median_size, median_tres=median_tres)

            # Check which bad pixels are new. Also, compare the
            # analytically determined noise with the empirically measured
            # noise.
            pxdq_new = np.sum(temp[pxdq_cut < 0.5])
            log.info(
                f"         Iteration {k + 1:.0f}: {pxdq_new:.0f} new bad pixels, "
                f"sdev of norm noise = {np.std(fmas_data[pmas]):.3f}"
            )

            # If no new bad pixels were identified, terminate the
            # iteration.
            if pxdq_new == 0:
                break

            # If new bad pixels were identified, add them to the bad pixel
            # map.
            pxdq_cut = ((pxdq_cut > 0.5) | (temp > 0.5)).astype("bool")

        # Put the modified frames back into the data cube.
        data[j, :idx_x, :idx_y] = fourier_corr(data_orig, pxdq_cut, fmas)
        pxdq[j, :idx_x, :idx_y] = pxdq_cut

    return data, pxdq
