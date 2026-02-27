import logging
import warnings

import numpy as np
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.utils.exceptions import AstropyUserWarning
from photutils.background import Background2D, MedianBackground
from scipy.optimize import curve_fit

log = logging.getLogger(__name__)

__all__ = ["clip_to_background", "background_level"]


def clip_to_background(
    image,
    mask,
    sigma_lower=3.0,
    sigma_upper=2.0,
    fit_histogram=False,
    lower_half_only=False,
    verbose=False,
):
    """
    Flag signal and bad pixels in the image mask.

    Given an image, estimate the background level and a sigma value for the
    mean background.

    The center and sigma may be calculated from a simple sigma-clipped
    median and standard deviation, or may be refined by fitting a Gaussian
    to a histogram of the values.  In that case, the center is the
    Gaussian mean and the sigma is the Gaussian width.

    Pixels above the ``center + sigma_upper * sigma`` are assumed to
    have signal; those below this level are assumed to be background pixels.

    Pixels less than ``center - sigma_lower * sigma`` are also excluded as bad values.

    The input mask is updated in place.

    Parameters
    ----------
    image : ndarray of float
        2D image containing signal and background values.
    mask : ndarray of bool
        2D input mask to be updated. True indicates background
        pixels to be used. Regions containing signal or outliers
        will be set to `False`.
    sigma_lower : float, optional
        The number of standard deviations to use as the lower bound
        for the clipping limit. Values below this limit are marked
        `False` in the mask.
    sigma_upper : float, optional
        The number of standard deviations to use as the upper bound
        for the clipping limit. Values above this limit are marked
        `False` in the mask.
    fit_histogram :  bool, optional
        If set, the center value and standard deviation used with
        ``sigma_lower`` and ``sigma_upper`` for clipping outliers is derived
        from a Gaussian fit to a histogram of values. Otherwise, the
        center and standard deviation are derived from a simple iterative
        sigma clipping.
    lower_half_only : bool, optional
        If set, the data used to compute the center and standard deviation
        for clipping is the lower half of the distribution only. Values
        below the median are mirrored around the median value to simulate
        a symmetric distribution.  This is intended to account for
        asymmetrical value distributions, with long tails in the upper
        half of the distribution, due to diffuse emission, for example.
    verbose : bool, optional
        If set, DEBUG level messages are issued with details on the
        computed statistics.
    """
    # Use float64 for stats computations
    image = image.astype(np.float64)

    # Sigma limit for basic stats
    sigma_limit = 3.0

    # Check mask for any valid data before proceeding
    if not np.any(mask):
        return

    # Initial iterative sigma clip
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning, message=".* slice")
        mean, median, sigma = sigma_clipped_stats(image, mask=~mask, sigma=sigma_limit)
    if fit_histogram:
        center = mean
    else:
        center = median
    if verbose:
        log.debug("From initial sigma clip:")
        log.debug(f"    center: {center:.5g}")
        log.debug(f"    sigma: {sigma:.5g}")

    # If desired, use only the lower half of the data distribution
    if lower_half_only:
        lower_half_idx = mask & (image <= center)
        data_for_stats = (
            np.concatenate(((image[lower_half_idx] - center), (center - image[lower_half_idx])))
            + center
        )

        # Redo stats on lower half of distribution
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
            mean, median, sigma = sigma_clipped_stats(data_for_stats, sigma=sigma_limit)
        if fit_histogram:
            center = mean
        else:
            center = median
        if verbose:
            log.debug("From lower half distribution:")
            log.debug(f"    center: {center:.5g}")
            log.debug(f"    sigma: {sigma:.5g}")
    else:
        data_for_stats = image[mask]

    # Refine sigma and center from a fit to a histogram, if desired
    if fit_histogram:
        try:
            hist, edges = np.histogram(
                data_for_stats, bins=2000, range=(center - 4.0 * sigma, center + 4.0 * sigma)
            )
        except ValueError:
            log.error("Histogram failed; using clip center and sigma.")
            hist, edges = None, None

        param_opt = None
        if hist is not None:
            values = (edges[1:] + edges[0:-1]) / 2.0
            ind = np.argmax(hist)
            mode_estimate = values[ind]

            # Fit a Gaussian profile to the histogram
            def gaussian(x, g_amp, g_mean, g_sigma):
                return g_amp * np.exp(-0.5 * ((x - g_mean) / g_sigma) ** 2)

            param_start = (hist[ind], mode_estimate, sigma)
            bounds = [(0, values[0], 0), (np.inf, values[-1], values[-1] - values[0])]
            try:
                param_opt, _ = curve_fit(gaussian, values, hist, p0=param_start, bounds=bounds)
            except RuntimeError:
                log.error("Gaussian fit failed; using clip center and sigma.")
                param_opt = None

            if verbose:
                log.debug("From histogram:")
                log.debug(f"    mode estimate: {mode_estimate:.5g}")
                log.debug(f"    range of values in histogram: {values[0]:.5g} to {values[-1]:.5g}")

        if verbose:
            log.debug("Gaussian fit results:")
        if param_opt is None:
            if verbose:
                log.debug("    (fit failed)")
        else:
            if verbose:
                log.debug(f"    peak: {param_opt[0]:.5g}")
                log.debug(f"    center: {param_opt[1]:.5g}")
                log.debug(f"    sigma: {param_opt[2]:.5g}")
            center = param_opt[1]
            sigma = param_opt[2]

    # Set limits from center and sigma
    background_lower_limit = center - sigma_lower * sigma
    background_upper_limit = center + sigma_upper * sigma
    if verbose:
        log.debug(f"Mask limits: {background_lower_limit:.5g} to {background_upper_limit:.5g}")

    # Clip bad values
    bad_values = image < background_lower_limit
    mask[bad_values] = False

    # Clip signal (> N sigma)
    signal = image > background_upper_limit
    mask[signal] = False


def background_level(image, mask, background_method="median", background_box_size=None):
    """
    Fit a low-resolution background level.

    Parameters
    ----------
    image : ndarray of float
        The 2D image containing the background to fit.
    mask : ndarray of bool
        The mask that indicates which pixels are to be used in fitting.
        True indicates a background pixel.
    background_method : {'median', 'model', None}, optional
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is fit with a low-resolution model via
        `~photutils.background.Background2D`.  If None, the background
        value is 0.0.
    background_box_size : tuple of int, optional
        Box size for the data grid used by `~photutils.background.Background2D` when
        ``background_method = 'model'``. For best results, use a box size
        that evenly divides the input image shape. Defaults to 32x32 if
        not provided.

    Returns
    -------
    background : float or ndarray of float
        The background level: a single value, if ``background_method``
        is 'median' or None, or an array matching the input image size
        if ``background_method`` is 'model'.
    """
    if background_method is None:
        background = 0.0

    else:
        # Sigma limit for basic stats
        sigma_limit = 3.0

        # Flag more signal in the background subtracted image,
        # with sigma set by the lower half of the distribution only
        clip_to_background(
            image, mask, sigma_lower=sigma_limit, sigma_upper=sigma_limit, lower_half_only=True
        )

        if background_method == "model":
            sigma_clip_for_bkg = SigmaClip(sigma=sigma_limit, maxiters=5)
            bkg_estimator = MedianBackground()

            if background_box_size is None:
                # use 32 x 32 if possible, otherwise take next largest box
                # size that evenly divides the image (minimum 1)
                background_box_size = []
                recommended = np.arange(1, 33)
                for i_size in image.shape:
                    divides_evenly = i_size % recommended == 0
                    background_box_size.append(int(recommended[divides_evenly][-1]))
                log.debug(f"Using box size {background_box_size}")

            box_division_remainder = (
                image.shape[0] % background_box_size[0],
                image.shape[1] % background_box_size[1],
            )
            if not np.allclose(box_division_remainder, 0):
                log.warning(
                    f"Background box size {background_box_size} "
                    f"does not divide evenly into the image "
                    f"shape {image.shape}."
                )

            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
                    bkg = Background2D(
                        image,
                        box_size=background_box_size,
                        filter_size=(5, 5),
                        mask=~mask,
                        sigma_clip=sigma_clip_for_bkg,
                        bkg_estimator=bkg_estimator,
                    )
                background = bkg.background
            except ValueError:
                log.error("Background fit failed, using median value.")
                background = np.nanmedian(image[mask])
        else:
            background = np.nanmedian(image[mask])
    return background
