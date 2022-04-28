"""
Module to detect sources using image segmentation.
"""

import logging
import warnings

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip
from astropy.utils import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from photutils.background import Background2D, MedianBackground
from photutils.utils.exceptions import NoDetectionsWarning
from photutils.segmentation import detect_sources, deblend_sources

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class JWSTBackground:
    """
    Class to estimate a 2D background and background RMS noise in an
    image.

    Parameters
    ----------
    data : 2D `~numpy.ndarray`
        The input 2D array.

    box_size : int or array_like (int)
        The box size along each axis.  If ``box_size`` is a scalar then
        a square box of size ``box_size`` will be used.  If ``box_size``
        has two elements, they should be in ``(ny, nx)`` order.

    coverage_mask : array_like (bool), optional
        A boolean mask, with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked data are excluded from calculations. ``coverage_mask``
        should be `True` where there is no coverage (i.e., no data) for
        a given pixel (e.g., blank areas in a mosaic image). It should
        not be used for bad pixels.

    Attributes
    ----------
    background : 2D `~numpy.ndimage`
        The estimated 2D background image.

    background_rms : 2D `~numpy.ndimage`
        The estimated 2D background RMS image.
    """

    def __init__(self, data, box_size=100, coverage_mask=None):
        self.data = data
        self.box_size = np.asarray(box_size).astype(int)  # must be integer
        self.coverage_mask = coverage_mask

    @lazyproperty
    def _background2d(self):
        """
        Estimate the 2D background and background RMS noise in an image.

        Returns
        -------
        background : `photutils.background.Background2D`
            A Background2D object containing the 2D background and
            background RMS noise estimates.
        """
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        filter_size = (3, 3)

        try:
            bkg = Background2D(self.data, self.box_size,
                               filter_size=filter_size,
                               coverage_mask=self.coverage_mask,
                               sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator)
        except ValueError:
            # use the entire unmasked array
            bkg = Background2D(self.data, self.data.shape,
                               filter_size=filter_size,
                               coverage_mask=self.coverage_mask,
                               sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator,
                               exclude_percentile=100.)
            log.info('Background could not be estimated in meshes. '
                     'Using the entire unmasked array for background '
                     f'estimation: bkg_boxsize={self.data.shape}.')

        return bkg

    @lazyproperty
    def background(self):
        """
        The 2D background image.
        """
        return self._background2d.background

    @lazyproperty
    def background_rms(self):
        """
        The 2D background RMS image.
        """
        return self._background2d.background_rms


def make_kernel(kernel_fwhm):
    """
    Make a 2D Gaussian smoothing kernel that is used to filter the image
    before thresholding.

    Filtering the image will smooth the noise and maximize detectability
    of objects with a shape similar to the kernel.

    The kernel must have odd sizes in both X and Y, be centered in the
    central pixel, and normalized to sum to 1.

    Parameters
    ----------
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.

    Returns
    -------
    kernel : `astropy.convolution.Kernel2D`
        The output smoothing kernel, normalized such that it sums to 1.
    """
    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma)
    kernel.normalize(mode='integral')
    return kernel


def convolve_data(data, kernel_fwhm, mask=None):
    """
    Convolve the data with a Gaussian2D kernel.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The 2D array to convolve.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.

    mask : array_like, bool, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
    """
    kernel = make_kernel(kernel_fwhm)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)
        return convolve(data, kernel, mask=mask, normalize_kernel=True)


class JWSTSourceFinder:
    """
    Class to detect sources, including deblending, using image
    segmentation.

    Parameters
    ----------
    threshold : float
        The data value to be used as the per-pixel detection threshold.

    npixels : int
        The number of connected pixels, each greater than the threshold,
        that an object must have to be detected. ``npixels`` must be a
        positive integer.

    deblend : bool, optional
        Whether to deblend overlapping sources. Source deblending
        requires scikit-image.
    """

    def __init__(self, threshold, npixels, deblend=False):
        self.threshold = threshold
        self.npixels = npixels
        self.deblend = deblend
        self.connectivity = 8
        self.nlevels = 32
        self.contrast = 0.001
        self.mode = 'exponential'

    def __call__(self, convolved_data, mask=None):
        """
        Parameters
        ----------
        convolved_data : 2D `numpy.ndarray`
            The 2D convolved array from which to detect sources.

        mask : array_like, bool, optional
            A boolean mask with the same shape as ``convolved_data``,
            where a `True` value indicates the corresponding element
            of ``convolved_data`` is masked. Masked pixels will not be
            included in any source.

        Returns
        -------
        segment_image : `~photutils.segmentation.SegmentationImage` or `None`
            A 2D segmentation image, with the same shape as the input data,
            where sources are marked by different positive integer values. A
            value of zero is reserved for the background. If no sources are
            found then `None` is returned.
        """
        if mask is not None:
            if mask.all():
                log.error('There are no valid pixels in the image to detect '
                          'sources.')
                return None

        with warnings.catch_warnings():
            # suppress NoDetectionsWarning from photutils
            warnings.filterwarnings('ignore', category=NoDetectionsWarning)

            segment_img = detect_sources(convolved_data, self.threshold,
                                         self.npixels, mask=mask,
                                         connectivity=self.connectivity)
            if segment_img is None:
                log.warning('No sources were found. Source catalog will not '
                            'be created.')
                return None

            # source deblending requires scikit-image
            if self.deblend:
                segment_img = deblend_sources(convolved_data, segment_img,
                                              npixels=self.npixels,
                                              nlevels=self.nlevels,
                                              contrast=self.contrast,
                                              mode=self.mode,
                                              connectivity=self.connectivity,
                                              relabel=True)

        log.info(f'Detected {segment_img.nlabels} sources')
        return segment_img
