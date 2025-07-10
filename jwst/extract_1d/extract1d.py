import warnings

import numpy as np
from astropy import convolution

__all__ = ["extract1d"]


def build_coef_matrix(image, profiles_2d=None, profile_bg=None, weights=None, order=0):
    """
    Build matrices and vectors to enable least-squares fits.

    Parameters
    ----------
    image : ndarray
        2D array, transposed if necessary so that the dispersion direction
        is the second index.
    profiles_2d : list of ndarray or None, optional
        These 2D arrays contain the weights for the extraction.  These arrays
        should be the same shape as image, with one array for each object
        to extract.  If set to None, the coefficient matrix values default
        to unity.
    profile_bg : ndarray of bool or None, optional
        2D array of the same shape as image, with nonzero elements where the
        background is to be estimated.  If not specified, no additional
        pixels are included for background calculations.
    weights : ndarray or None, optional
        2D array of (float) weights for the extraction.  If using inverse
        variance weighting, these should be the square root of the inverse
        variance.  If not supplied, unit weights will be used.
    order : int, optional
        Polynomial order for fitting to each column of background.
        Default 0 (uniform background).

    Returns
    -------
    matrix : ndarray, 3-D, float64
        Design matrix for each pixel, shape (npixels, npar, npar)
    vec : ndarray, 2-D, float64
        Target vectors for the design matrix, shape (npixels, npar)
    coefmatrix : ndarray, 3-D, float64
        Matrix of coefficients for each parameter for each pixel,
        used to reconstruct the model fit.  Shape (npixels, npixels_y, npar)
    """
    if profiles_2d is None:
        profiles_2d = []

    # Independent variable values for the polynomial fit.
    y = np.linspace(-1, 1, image.shape[0])

    # Build the matrix of terms that multiply the coefficients.
    # Polynomial terms first, then source terms if those arrays
    # are supplied.
    coefmatrix = np.ones((image.shape[1], image.shape[0], order + 1 + len(profiles_2d)))
    for i in range(1, order + 1):
        coefmatrix[..., i] = coefmatrix[..., i - 1] * y[np.newaxis, :]

    # Here are the source terms.
    for i in range(len(profiles_2d)):
        coefmatrix[..., i + order + 1] = profiles_2d[i].T

    # Construct a boolean array for the pixels that are nonzero
    # in any of our profiles.
    pixels_used = np.zeros(image.T.shape, dtype=bool)
    for profile_2d in profiles_2d:
        pixels_used = pixels_used | (profile_2d.T != 0)
    if profile_bg is not None:
        pixels_used = pixels_used | (profile_bg.T != 0)
    pixels_used = pixels_used & np.isfinite(image.T)

    # Target vector and coefficient vector for the least squares fit.
    targetvector = image.T.copy()

    # We don't want to be ruined by NaNs in regions we are not fitting anyway.
    targetvector[~pixels_used] = 0
    coefmatrix_masked = coefmatrix * pixels_used[:, :, np.newaxis]

    # Weighting goes here.  If we are using inverse variance weighting,
    # weight here is the square root of the inverse variance.
    if weights is not None:
        coefmatrix_masked *= weights.T[:, :, np.newaxis]
        targetvector *= weights.T

    # Products of the coefficient matrices suitable for passing to
    # linalg.solve.  These are matrices of size (npixels, npar, npar)
    # and (npixels, npar).
    matrix = np.einsum("lji,ljk->lik", coefmatrix_masked, coefmatrix_masked)  # codespell:ignore
    vec = np.einsum("lji,lj->li", coefmatrix_masked, targetvector)

    return matrix, vec, coefmatrix, coefmatrix_masked


def _fit_background_for_box_extraction(
    image,
    profiles_2d,
    variance_rn,
    variance_phnoise,
    variance_flat,
    profile_bg,
    weights,
    bg_smooth_length,
    bkg_fit_type,
    bkg_order,
):
    """
    Fit a background level for box extraction.

    Returns
    -------
    bkg_2d, var_bkg_rn, var_bkg_phnoise, var_bkg_flat : tuple of ndarray
        Background and associated variances.
    """
    # Start by copying the input image
    input_background = image.copy()

    # Smooth the image, if desired, for computing a background.
    # Astropy's convolve routine will replace NaNs.  The convolution
    # is done along the dispersion direction.
    if bg_smooth_length > 1:
        if not bg_smooth_length % 2 == 1:
            raise ValueError("bg_smooth_length should be an odd integer >= 1.")
        kernel = np.ones((1, bg_smooth_length)) / bg_smooth_length
        input_background = convolution.convolve(input_background, kernel, boundary="extend")

    if bkg_fit_type == "median":
        input_background[profile_bg == 0] = np.nan
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")
            bkg_1d = np.nanmedian(input_background, axis=0)
        bkg_2d = bkg_1d[np.newaxis, :]

        # Putting an uncertainty on the median is a bit harder.
        # It is typically about 1.2 times the uncertainty on the mean.
        # The details depend on the number of pixels averaged...
        # bkg_npix is the total weight that we will need to apply
        # to the background value when removing it from the 2D image.
        # pixwgt normalizes the total weights of all pixels used to 1.
        bkg_npix = np.sum(profiles_2d[0], axis=0)
        variance = variance_rn + variance_phnoise + variance_flat
        wgt = np.isfinite(image) * np.isfinite(variance) * (profile_bg != 0)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value")
            pixwgt = wgt / np.sum(wgt, axis=0)[np.newaxis, :]

        var_bkg_rn = 1.2**2 * bkg_npix**2 * np.array([np.nansum(variance_rn * pixwgt**2, axis=0)])
        var_bkg_phnoise = (
            1.2**2 * bkg_npix**2 * np.array([np.nansum(variance_phnoise * pixwgt**2, axis=0)])
        )
        var_bkg_flat = (
            1.2**2 * bkg_npix**2 * np.array([np.nansum(variance_flat * pixwgt**2, axis=0)])
        )

    elif bkg_fit_type == "poly" and bkg_order >= 0:
        if not bkg_order == int(bkg_order):
            raise ValueError("If bkg_fit_type is 'poly', bkg_order must be an integer >= 0.")

        # Build the matrices to fit a polynomial column-by-column.
        result = build_coef_matrix(
            input_background, profile_bg=profile_bg, weights=weights, order=bkg_order
        )
        matrix, vec, coefmatrix, coefmatrix_masked = result

        # Don't try to solve singular matrices.  Background will be
        # zero in these cases.  Could make them NaN if you want.
        ok = np.linalg.cond(matrix) < 1e10
        cov_bg_coefs = np.zeros(matrix.shape)
        cov_bg_coefs[ok] = np.linalg.inv(matrix[ok])

        # These are the pixel-dependent weights to compute our coefficients.
        # We will use them to propagate errors.
        pixwgt = weights.T[:, :, np.newaxis] * np.einsum(
            "ijk,ilj->ilk", cov_bg_coefs, coefmatrix_masked
        )
        bkg_mat = np.sum(np.swapaxes(coefmatrix, 0, 1) * profiles_2d[0][:, :, np.newaxis], axis=0)

        # Sum of all the contributions to the background at the pixels
        # where we will do the extraction.  Used to propagate errors.
        pixwgt_tot = np.sum(bkg_mat[:, np.newaxis, :] * pixwgt, axis=-1)

        var_bkg_rn = np.array([np.nansum(variance_rn.T * pixwgt_tot**2, axis=1)])
        var_bkg_phnoise = np.array([np.nansum(variance_phnoise.T * pixwgt_tot**2, axis=1)])
        var_bkg_flat = np.array([np.nansum(variance_flat.T * pixwgt_tot**2, axis=1)])

        coefs = np.einsum("ijk,ij->ik", cov_bg_coefs, vec)

        # Reconstruct the 2D background.
        bkg_2d = np.sum(coefs[:, np.newaxis, :] * coefmatrix, axis=-1).T

    else:
        raise ValueError(
            "bkg_fit_type should be 'median' or 'poly'. "
            "If 'poly', bkg_order must be an integer >= 0."
        )

    return bkg_2d, var_bkg_rn, var_bkg_phnoise, var_bkg_flat


def _box_extract(
    image,
    profiles_2d,
    variance_rn,
    variance_phnoise,
    variance_flat,
    bkg_2d,
    var_bkg_rn,
    var_bkg_phnoise,
    var_bkg_flat,
    model,
):
    """
    Perform box extraction.

    Returns
    -------
    tuple of ndarray
        Extracted spectrum and associated data.
    """
    # This only makes sense with a single profile, i.e., pulling out
    # a single spectrum.
    nobjects = 1
    profile_2d = profiles_2d[0].copy()
    image_masked = image.copy()

    # Mask NaNs and infs for the extraction
    profile_2d[~np.isfinite(image_masked)] = 0
    image_masked[profile_2d == 0] = 0

    # Return array of shape (1, npixels) for generality
    fluxes = np.array([np.sum(image_masked * profile_2d, axis=0)])

    # Number of contributing pixels at each wavelength.
    npixels = np.array([np.sum(profile_2d, axis=0)])

    # Add average flux over the aperture to the model, so that
    # a sum over the cross-dispersion direction reproduces the summed flux
    valid = npixels[0] > 0
    model[:, valid] += (fluxes[0][valid] / npixels[0][valid]) * profile_2d[:, valid]

    # Compute the variance on the sum, same shape as f.
    # Need to decompose this into read noise, photon noise, and flat noise.

    # Ignore overflow warnings and invalid values - some variance values may be
    # very small or very large.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "overflow encountered", RuntimeWarning)
        warnings.filterwarnings("ignore", "invalid value", RuntimeWarning)
        var_rn = np.array([np.nansum(variance_rn * profile_2d**2, axis=0)])
        var_phnoise = np.array([np.nansum(variance_phnoise * profile_2d**2, axis=0)])
        var_flat = np.array([np.nansum(variance_flat * profile_2d**2, axis=0)])

    if bkg_2d is not None:
        var_rn += var_bkg_rn
        var_phnoise += var_bkg_phnoise
        var_flat += var_bkg_flat
        bkg = np.array([np.sum(bkg_2d * profile_2d, axis=0)])
        bkg[~np.isfinite(bkg)] = 0.0
    else:
        bkg = np.zeros((nobjects, image.shape[1]))

    return (
        fluxes,
        var_rn,
        var_phnoise,
        var_flat,
        bkg,
        var_bkg_rn,
        var_bkg_phnoise,
        var_bkg_flat,
        npixels,
        model,
    )


def _optimal_extract(
    image,
    profiles_2d,
    variance_rn,
    variance_phnoise,
    variance_flat,
    weights,
    profile_bg,
    fit_bkg,
    bkg_order,
    bkg_2d,
    var_bkg_rn,
    var_bkg_phnoise,
    var_bkg_flat,
    model,
):
    """
    Perform optimal extraction.

    Returns
    -------
    tuple of ndarray
        Extracted spectrum and associated data.
    """
    # Background fitting needs to be done simultaneously with the
    # fitting of the spectra in this case.  If we are not fitting a
    # background, pass -1 for the order of the polynomial correction.
    if fit_bkg:
        order = bkg_order
        if order != int(order) or order < 0:
            raise ValueError("For optimal extraction, bkg_order must be an integer >= 0.")
    else:
        order = -1

    result = build_coef_matrix(
        image, profiles_2d=profiles_2d, weights=weights, profile_bg=profile_bg, order=order
    )
    matrix, vec, coefmatrix, coefmatrix_masked = result

    # Don't try to solve equations with singular matrices.
    # Fluxes will be zero in these cases.  We will make them NaN later.
    ok = np.linalg.cond(matrix) < 1e10

    # These are the covariance matrices for all parameters if inverse
    # variance weights are passed to the build_coef_matrix above.  For
    # generality, variances are actually computed using the weights on
    # each pixel and the associated variance in the input image.
    covariances = np.zeros(matrix.shape)
    covariances[ok] = np.linalg.inv(matrix[ok])

    # These are the pixel-dependent weights to compute our coefficients.
    # We will use them to propagate errors.
    pixwgt = weights.T[:, :, np.newaxis] * np.einsum("ijk,ilj->ilk", covariances, coefmatrix_masked)

    # Don't use NaN pixels in the sum.  These will already be zero in
    # pixwgt.  coefs are the best-fit coefficients of the source and
    # background components.
    coefs = np.nansum(pixwgt * image.T[:, :, np.newaxis], axis=1)

    # Effective number of contributing pixels at each wavelength for each source.
    nobjects = len(profiles_2d)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value")
        wgt_src_pix = [
            profiles_2d[i] * (weights > 0) / np.sum(profiles_2d[i] ** 2, axis=0)
            for i in range(nobjects)
        ]
    npixels = np.sum(wgt_src_pix, axis=1)

    if order > -1:
        bkg_2d = np.sum(coefs[:, np.newaxis, : order + 1] * coefmatrix[..., : order + 1], axis=-1).T

    # Variances for each object (discard variances for background here)
    var_rn = np.nansum(pixwgt[..., -nobjects:] ** 2 * variance_rn.T[:, :, np.newaxis], axis=1).T
    var_phnoise = np.nansum(
        pixwgt[..., -nobjects:] ** 2 * variance_phnoise.T[:, :, np.newaxis], axis=1
    ).T
    var_flat = np.nansum(pixwgt[..., -nobjects:] ** 2 * variance_flat.T[:, :, np.newaxis], axis=1).T

    # Computing a background contribution to the noise is harder in a joint fit.
    # Here, I am computing the weighting coefficients I would have without a background.
    # I then compute those variances, and subtract them from the actual variances.
    if order > -1:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value")
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero")

            wgt_nobkg = [
                profiles_2d[i] * weights / np.sum(profiles_2d[i] ** 2 * weights, axis=0)
                for i in range(nobjects)
            ]

            bkg = np.array([np.sum(wgt_nobkg[i] * bkg_2d, axis=0) for i in range(nobjects)])

            # Avoid overflow in squaring weights by multiplying by variance first
            var_bkg_rn = np.array(
                [
                    var_rn[i] - np.sum(variance_rn * wgt_nobkg[i] * wgt_nobkg[i], axis=0)
                    for i in range(nobjects)
                ]
            )
            var_bkg_phnoise = np.array(
                [
                    var_phnoise[i] - np.sum(variance_phnoise * wgt_nobkg[i] * wgt_nobkg[i], axis=0)
                    for i in range(nobjects)
                ]
            )
            var_bkg_flat = np.array(
                [
                    var_flat[i] - np.sum(variance_flat * wgt_nobkg[i] * wgt_nobkg[i], axis=0)
                    for i in range(nobjects)
                ]
            )

        # Make sure background values are finite
        bkg[~np.isfinite(bkg)] = 0.0

        # We did our best to estimate the background contribution to the variance
        # in this case.  Don't let it go negative.
        var_bkg_rn *= var_bkg_rn > 0
        var_bkg_phnoise *= var_bkg_phnoise > 0
        var_bkg_flat *= var_bkg_flat > 0
    else:
        bkg = np.zeros((nobjects, image.shape[1]))

    # Reshape to (nobjects, npixels)
    fluxes = coefs[:, -nobjects:].T
    model += np.sum(coefs[:, np.newaxis, :] * coefmatrix, axis=-1).T

    return (
        fluxes,
        var_rn,
        var_phnoise,
        var_flat,
        bkg,
        var_bkg_rn,
        var_bkg_phnoise,
        var_bkg_flat,
        npixels,
        model,
    )


def extract1d(
    image,
    profiles_2d,
    variance_rn,
    variance_phnoise,
    variance_flat,
    weights=None,
    profile_bg=None,
    extraction_type="box",
    bg_smooth_length=0,
    fit_bkg=False,
    bkg_fit_type="poly",
    bkg_order=0,
):
    """
    Extract the spectrum, optionally subtracting background.

    Parameters
    ----------
    image : ndarray
        2D array, transposed if necessary so that the dispersion direction
        is the second index.
    profiles_2d : list of ndarray
        These 2D arrays contain the weights for the extraction.  A box
        extraction will add up the flux multiplied by these weights; an
        optimal extraction will fit an amplitude to the weight map at each
        column in the dispersion direction.  These arrays should be the
        same shape as image, with one array for each object to extract.
        Box extraction only works if exactly one profile is supplied
        (i.e. this is a one-element list).
    variance_rn : ndarray
        2D read noise component of the variance.
    variance_phnoise : ndarray
        2D photon noise component of the variance.
    variance_flat : ndarray
        2D flat component of the variance.
    weights : ndarray or None, optional
        2D weights for the individual pixels in fitting a profile.  If None
        (default), use uniform weights (ones for all valid pixels).
    profile_bg : ndarray or None, optional
        2D array of the same shape as image, with nonzero elements where the
        background is to be estimated.
    extraction_type : {"box", "optimal"}, optional
        Type of spectral extraction.
    bg_smooth_length : int, optional
        Smoothing length for box smoothing of the background along the
        dispersion direction.  Should be odd, >=1.
    fit_bkg : bool, optional
        Fit a background?  Default False
    bkg_fit_type : str, optional
        Type of fitting to apply to background values in each column (or
        row, if the dispersion is vertical).
    bkg_order : int, optional
        Polynomial order for fitting to each column of background.  A value
        of 0 means that a simple average of the background regions, column
        by column, will be used.
        This argument must be positive or zero, and it is only used if
        background regions have been specified and if `bkg_fit` is `poly`.

    Returns
    -------
    fluxes : ndarray of float64
        The extracted spectrum/spectra.  Units are currently arbitrary.
        The first dimension is the same as the length of profiles_2d.
    var_rn : ndarray of float64
        The variances of the extracted spectrum/spectra due to read noise.
        Units are the same as flux^2, shape is the same as flux.
    var_phnoise : ndarray of float64
        The variances of the extracted spectrum/spectra due to photon noise.
        Units are the same as flux^2, shape is the same as flux.
    var_flat : ndarray of float64
        The variances of the extracted spectrum/spectra due to flatfield
        uncertainty. Units are the same as flux^2, shape is the same as flux.
    bkg : ndarray of float64
        Background level that would be obtained for each source if performing
        a 1-D extraction on the 2D background.
    var_bkg_rn : ndarray of float64
        As above, for read noise.  Nonzero because read noise adds an error term
        to the derived background level.
    var_bkg_phnoise : ndarray of float64
        The variances of the extracted spectrum/spectra due to background photon
        noise. Units are the same as flux^2, shape is the same as flux.  This
        background contribution is already included in var_phnoise.
    var_bkg_flat : ndarray of float64
        As above, for the flatfield.
    npixels : ndarray of int64
        Number of pixels that contribute to the flux measurement for each source
    model : ndarray of float64
        The model of the scene, the same shape as the input image (and
        hopefully also similar in value).
    """
    nobjects = len(profiles_2d)  # hopefully at least one!
    model = np.zeros(image.shape)

    # If weights are not supplied, equally weight all valid pixels.
    variance = variance_rn + variance_phnoise + variance_flat
    if weights is None:
        weights = np.isfinite(variance) * np.isfinite(image)

    # This is the case of a background fit independent of a flux fit.
    # This is done only if we have a background region to fit, the
    # boolean variable to fit is set, and we are using box extraction.
    # Inverse variance weights should be used with care, as they have the
    # potential to introduce biases.
    if profile_bg is not None and fit_bkg and extraction_type == "box":
        (bkg_2d, var_bkg_rn, var_bkg_phnoise, var_bkg_flat) = _fit_background_for_box_extraction(
            image,
            profiles_2d,
            variance_rn,
            variance_phnoise,
            variance_flat,
            profile_bg,
            weights,
            bg_smooth_length,
            bkg_fit_type,
            bkg_order,
        )

        model += bkg_2d
        image_sub = image - bkg_2d
    else:
        image_sub = image.copy()

        # Set background image to None
        bkg_2d = None

        # Set background uncertainties to zero.
        var_bkg_rn = np.zeros((nobjects, image.shape[1]))
        var_bkg_phnoise = np.zeros((nobjects, image.shape[1]))
        var_bkg_flat = np.zeros((nobjects, image.shape[1]))

    # This is the case of box extraction.
    if extraction_type == "box" and len(profiles_2d) == 1:
        (
            fluxes,
            var_rn,
            var_phnoise,
            var_flat,
            bkg,
            var_bkg_rn,
            var_bkg_phnoise,
            var_bkg_flat,
            npixels,
            model,
        ) = _box_extract(
            image_sub,
            profiles_2d,
            variance_rn,
            variance_phnoise,
            variance_flat,
            bkg_2d,
            var_bkg_rn,
            var_bkg_phnoise,
            var_bkg_flat,
            model,
        )

    # This is optimal extraction (well, profile-based extraction).  It
    # is "optimal extraction" if the weights are the inverse variances,
    # but you must be careful about biases in this case.
    elif extraction_type == "optimal":
        (
            fluxes,
            var_rn,
            var_phnoise,
            var_flat,
            bkg,
            var_bkg_rn,
            var_bkg_phnoise,
            var_bkg_flat,
            npixels,
            model,
        ) = _optimal_extract(
            image_sub,
            profiles_2d,
            variance_rn,
            variance_phnoise,
            variance_flat,
            weights,
            profile_bg,
            fit_bkg,
            bkg_order,
            bkg_2d,
            var_bkg_rn,
            var_bkg_phnoise,
            var_bkg_flat,
            model,
        )

    else:
        raise ValueError(
            f"Extraction method {extraction_type} not supported with {nobjects} input profiles."
        )

    no_data = np.isclose(npixels, 0)
    fluxes[no_data] = np.nan

    return (
        fluxes,
        var_rn,
        var_phnoise,
        var_flat,
        bkg,
        var_bkg_rn,
        var_bkg_phnoise,
        var_bkg_flat,
        npixels,
        model,
    )
