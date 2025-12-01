import logging
import warnings

import gwcs
import numpy as np
from astropy.modeling.models import Identity, Scale, Shift
from astropy.stats import sigma_clipped_stats as scs
from astropy.utils.exceptions import AstropyUserWarning
from scipy.interpolate import make_lsq_spline
from scipy.signal import find_peaks

from jwst.assign_wcs.nirspec import nrs_ifu_wcs

__all__ = ["bspline_fit", "fit_2d_spline_profile"]

log = logging.getLogger(__name__)


def bspline_fit(
    xvec,
    yvec,
    nbkpts=50,
    wrapsig_low=3.0,
    wrapsig_high=3.0,
    wrapiter=3,
    spaceratio=1.6,
    verbose=False,
):
    """
    Fit a univariate basis spline to a vector with iterative rejection.

    Parameters
    ----------
    xvec : ndarray
        1D array of x values for the vector.
    yvec : ndarray
        1D array of y values to fit.
    nbkpts : int, optional
        Number of spline breakpoints (knots).
    wrapsig_low : float, optional
        Low sigma threshold for iterative fit.
    wrapsig_high : float, optional
        High sigma threshold for iterative fit.
    wrapiter : int, optional
        Number of iterations for the fit.
    spaceratio : float, optional
        Maximum spacing ratio to allow fitting to continue. If
        the tenth-largest spacing in the input ``xvec`` is larger
        than the knot spacing by this ratio, then return None instead
        of attempting to fit the data.
    verbose : bool, optional
        If True, a debug log message is printed for each fitting iteration.

    Returns
    -------
    spline_model : `~scipy.interpolate.BSpline` or None
        The spline model fit to the data.  If the fit failed, None is returned.
    """
    spline_model = None
    indx = np.isfinite(xvec) & np.isfinite(yvec)
    xvec_use = xvec[indx]
    yvec_use = yvec[indx]

    # Reject points with abnormally large spacing as they can
    # cause problems with breakpoints
    spacing = np.abs(np.diff(xvec_use, prepend=0))
    space_mean, _, space_rms = scs(spacing)
    if space_rms > 0:
        good = spacing < (space_mean + 5 * space_rms)
        xvec_use = xvec_use[good]
        yvec_use = yvec_use[good]

    # If the tenth-largest spacing was bigger than the knot spacing
    # (with some margin) then don't bspline
    # This factor of 1.6 was dialed based on inspection of the results
    # as sampling gets progressively worse for NIRSpec detectors
    knotspacing = (np.max(xvec_use) - np.min(xvec_use)) / nbkpts
    tenthspace = np.partition(spacing, -10)[-10]
    if tenthspace > (spaceratio * knotspacing):
        return spline_model

    # Number of points before iterative loop
    norig = len(xvec_use)

    last_spline = None
    for ii in range(0, wrapiter):
        knotspacing = (np.max(xvec_use) - np.min(xvec_use)) / nbkpts
        knotmin = np.min(xvec_use)
        knotmax = np.max(xvec_use) - knotspacing / 2.0
        knots = np.arange(knotmin, knotmax, knotspacing)

        # Add exterior knots for cubic spline
        degree = 3
        xb = xvec_use[0]
        xe = xvec_use[-1]
        knots = np.concatenate(([xb] * (degree + 1), knots, [xe] * (degree + 1)))

        spline_model = make_lsq_spline(xvec_use, yvec_use, knots, k=degree)
        ytemp = spline_model(xvec_use)
        if np.any(np.isnan(ytemp)):
            # Bad fit: break the loop and return the last fit spline (or None)
            spline_model = last_spline
            break

        # Good fit: keep it for the next loop.
        last_spline = spline_model

        # Reject significant outliers
        diff = yvec_use - ytemp
        rms = np.std(diff)
        rej = (diff < -wrapsig_low * rms) | (diff > wrapsig_high * rms)
        keep = ~rej
        if verbose:
            log.debug(f"Iter {ii} Rejected {np.sum(rej)} Kept {np.sum(keep)}")

        # If no points rejected or too few kept break out of the loop
        if np.sum(keep) < 0.8 * norig:
            # Too many rejected - keep the last model instead
            # todo - instead return None?
            spline_model = last_spline
            break
        elif np.sum(rej) == 0:
            break

        # Otherwise update the x and y vectors and try again
        xvec_use = xvec_use[keep]
        yvec_use = yvec_use[keep]

    return spline_model


def _get_weights_for_fit_scale(ratio, model_fit):
    """
    Get weights for scaling the spline model to the fit data.

    Sets weight to zero for outliers and invalid data points.

    Parameters
    ----------
    ratio : ndarray
        Ratio between the fit data and the model evaluated at the data points.
    model_fit : ndarray
        The spline model evaluated at the data points.

    Returns
    -------
    weights : ndarray
        Array matching the ``ratio`` shape, containing weights for each
        ratio data point.
    """
    # Weights start off proportional to flux of the model
    weights = model_fit.copy()

    # Weights are zero for any pixel where the data was NaN
    weights[~np.isfinite(ratio)] = 0

    # Weights are zero for any pixel where the data or model was negative
    weights[(ratio < 0) | (model_fit < 0)] = 0

    # Identify the 5 largest weight points
    order = np.argsort(weights)
    largest_5 = (weights >= weights[order[-5]]) & (np.isfinite(ratio))

    # Sigma-clipped mean and rms of these 5 ratios
    mean, _, rms = scs(ratio[largest_5])

    # Bad if over 2 sigma away
    bad = np.abs(mean - ratio) > (2 * rms)
    weights[bad] = 0

    # Normalize weights
    weights /= np.nansum(weights)

    return weights


def fit_2d_spline_profile(
    flux,
    alpha,
    lrange=50,
    col_index=None,
    require_ngood=8,
    splinebkpt=36,
    spaceratio=1.2,
    fit_scale=None,
):
    """
    Create a profile from spline fits to a single slit/slice image.

    Image must be oriented so that wavelengths are along x-axis. Each
    column is fit separately, with a window to include nearby data.

    Parameters
    ----------
    flux : ndarray
        Input 2D flux image to fit.
    alpha : ndarray
        Alpha coordinates for input flux.
    lrange : int, optional
        Local column range for data to include in the fit, to the
        left and right of each input column.
    col_index : iterable or None, optional
        Iterable or generator that produces column index values to fit.
        If provided, columns will be fit in the order specified.
        If not provided, columns will be fit left to right.
    require_ngood : int, optional
        Minimum number of data points required to attempt a fit in a column.
    splinebkpt : int, optional
        Number of spline breakpoints (knots).
    spaceratio : float, optional
        Maximum spacing ratio to allow fitting to continue. If
        the tenth-largest spacing in the input ``xvec`` is larger
        than the knot spacing by this ratio, then return None instead
        of attempting to fit the data.
    fit_scale : ndarray, optional
        Array of scale values to apply to the input flux before fitting.

    Returns
    -------
    splines : dict
        Keys are column index numbers, values are `~scipy.interpolate.BSpline`.
        If a spline model could not be fit, the column index number is not present.
    scales : dict
        Keys are column index numbers, values are floating point scales, to pair with
        the returned models. If a spline model could not be fit, the column index
        number is not present.
    """
    # Define a fallback spline model, initialize to None
    spline_model_save = None

    # Set up the column fitting order if not provided
    xsize = flux.shape[-1]
    if col_index is None:
        col_index = range(0, xsize, 1)

    # Scale the flux for fitting
    if fit_scale is not None:
        scaled_flux = flux / fit_scale
    else:
        scaled_flux = flux

    # Loop over columns in the slit/slice
    splines = {}
    scales = {}
    for i in col_index:
        col_flux = flux[:, i]
        ngood = np.sum(np.isfinite(col_flux))
        if ngood <= require_ngood:
            continue

        # Get local alpha and flux values for fitting
        lstart = np.max([i - lrange, 0])
        lstop = np.min([i + lrange, xsize])
        local_alpha = alpha[:, lstart:lstop]
        local_data = scaled_flux[:, lstart:lstop]

        # Trim to finite values
        finite_alpha = np.isfinite(local_alpha)
        local_alpha = local_alpha[finite_alpha]
        local_data = local_data[finite_alpha]

        # Sort by alpha
        idx = np.argsort(local_alpha)
        local_alpha = local_alpha[idx]
        local_data = local_data[idx]

        # Fit a bspline to the local data
        try:
            bspline = bspline_fit(
                local_alpha,
                local_data,
                nbkpts=splinebkpt,
                wrapsig_low=2.5,
                wrapsig_high=2.5,
                wrapiter=3,
                spaceratio=spaceratio,
                verbose=False,
            )

            # If this routine could not get a fit (returned None) use the saved fit
            if bspline is None:
                spline_model = spline_model_save
            else:
                spline_model = bspline
        except Exception as err:
            log.warning(f"Spline fit failed at column {i}: {str(err)}")
            spline_model = spline_model_save

        # Check for a good model
        if spline_model is None:
            continue

        # Store the spline model for the column
        splines[i] = spline_model
        spline_model_save = spline_model

        # Evaluate the bspline at the valid input locations to determine
        # a scale factor for the fit
        col_alpha = alpha[:, i]
        idx = np.isfinite(col_alpha)
        col_alpha = col_alpha[idx]
        col_flux = col_flux[idx]
        col_fit = spline_model(col_alpha)

        # Determine the normalization by the weighted mean ratio between model and data
        # Weights are based on the model so that we can reject outliers
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=RuntimeWarning)

            ratio = col_flux / col_fit

            # Weights start off proportional to flux
            weights = _get_weights_for_fit_scale(ratio, col_fit)
            wmeanratio = np.nansum(ratio * weights)

        # Store the scale factor for the column
        scales[i] = wmeanratio

    return splines, scales


def _reindex(xmin, xmax, scale=2.0):
    """
    Convert pixel positions on the old grid to oversampled positions.

    For example, with oversample=2, [0, 1, 2] goes to
    [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25]. With oversample=3,
    [0, 1, 2] goes to [-0.33, 0, 0.33, 0.67, 1, 1.33, 1.67, 2, 2.33].

    Parameters
    ----------
    xmin : int
        Minimum index.
    xmax : int
        Maximum index.
    scale : float, optional
        Oversample scaling factor.

    Returns
    -------
    new_x : ndarray of int
        Array of indices in the new grid.
    old_x : ndarray of float
        Array of coordinates in the old grid, corresponding to the new indices.
    """
    # Indices in the new array
    new_x = np.arange(xmin * scale, (xmax + 1) * scale, dtype=np.int32)

    # Indices in the old array, scaled for new pixel spacing
    # Also offset to center new coordinates on the old
    old_x = new_x / scale - (scale - 1) / (scale * 2)
    return new_x, old_x


def profile_image(shape, spline_models, spline_scales, slice_map, alpha):
    profile = np.full(shape, np.nan, dtype=np.float32)
    alpha_slice = np.full(shape, np.nan, dtype=np.float32)
    profile_slice = np.full(shape, np.nan, dtype=np.float32)
    for slnum in spline_models:
        splines = spline_models[slnum]
        scales = spline_scales[slnum]
        alpha_slice[:] = np.nan
        profile_slice[:] = np.nan

        indx = slice_map == slnum
        alpha_slice[indx] = alpha[indx]

        # loop over columns
        for i in range(shape[-1]):
            if i not in splines:
                continue

            # Spline fit for relevant data
            alpha_col = alpha_slice[:, i]
            valid_alpha = np.isfinite(alpha_col)
            spline_fit = scales[i] * splines[i](alpha_col[valid_alpha])
            profile_slice[:, i][valid_alpha] = spline_fit

        profile[indx] = profile_slice[indx]

    return profile


def linear_oversample(data, slice_map, oversample_factor, require_ngood):
    slice_numbers = np.unique(slice_map)
    slice_numbers = slice_numbers[np.isfinite(slice_numbers)]

    ysize, xsize = data.shape
    _, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))

    os_shape = (int(np.ceil(ysize * oversample_factor)), xsize)
    os_data = np.full(os_shape, np.nan, dtype=np.float32)
    edge_limit = int(oversample_factor)

    data_slice = np.full_like(data, np.nan)
    y_slice = np.full_like(data, np.nan)
    for slnum in slice_numbers:
        data_slice[:] = np.nan
        y_slice[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = slice_map == slnum
        data_slice[indx] = data[indx]
        y_slice[indx] = basey[indx]

        for ii in range(xsize):
            valid_data = np.isfinite(data_slice[:, ii])
            ngood = np.sum(valid_data)
            if ngood <= require_ngood:
                continue

            col_y = y_slice[:, ii]
            col_flux = data_slice[:, ii]

            valid_y = np.isfinite(col_y)
            newy, oldy = _reindex(
                int(col_y[valid_y].min()), int(col_y[valid_y].max()), scale=oversample_factor
            )

            valid_interp = valid_y & valid_data
            interpval = np.interp(oldy, col_y[valid_interp], col_flux[valid_interp])
            if edge_limit >= 1:
                interpval[0:edge_limit] = np.nan
                interpval[-edge_limit:] = np.nan

            # Check for NaNs in the input: they should be preserved in the output
            closest_pix = np.round(oldy).astype(int)
            is_nan = ~np.isfinite(col_flux[closest_pix])
            interpval[is_nan] = np.nan

            os_data[newy, ii] = interpval

    return os_data


def just_fit(
    flux_orig,
    alpha_orig,
    slmap,
    detector,
    thresh,
    slstart=0,
    slstop=30,
    lrange=50,
    require_ngood=8,
    splinebkpt=36,
    spaceratio=1.2,
):
    ysize, xsize = flux_orig.shape

    # Set up the column fitting order by detector
    if detector == "NRS1":
        # For NRS1, start on the left of detector since the tilt wrt pixels is greatest here
        col_index = range(0, xsize, 1)
    elif detector == "NRS2":
        # For NRS2, start on the right of detector since the tilt wrt pixels is greatest here
        col_index = range(xsize - 1, 0, -1)
    elif detector.startswith("MIR"):
        # For MIRI we need to start on the left and run to the middle,
        # and then on the right to the middle in order to have the middle
        # section not go too far beyond last good fit
        col_index = np.concatenate(
            [np.arange(0, xsize // 2 + 1), np.arange(xsize - 1, xsize // 2, -1)]
        )
    else:
        raise ValueError("Unknown detector")

    # Arrays to reset with NaNs for each slice
    data_slice = np.full_like(flux_orig, np.nan)
    alpha_slice = np.full_like(flux_orig, np.nan)

    spline_models = {}
    spline_scales = {}
    for slnum in range(slstart, slstop):
        log.info("Fitting slice %s", slnum)

        # Reset holding arrays to NaN
        data_slice[:] = np.nan
        alpha_slice[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = slmap == slnum
        data_slice[indx] = flux_orig[indx]
        alpha_slice[indx] = alpha_orig[indx]

        # A running sum in a given detector column (used for normalization)
        runsum = np.nansum(data_slice, axis=0)

        # Collapse the slice along Y to get max in each column
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=RuntimeWarning)
            collapse = np.nanmax(data_slice, axis=0)

        # Median column max across all columns
        medcmax = np.nanmedian(collapse)

        # Is medcmax over threshold?  If so do, bspline for this slice.
        dospline = False
        if medcmax > thresh[slnum]:
            dospline = True

        if dospline:
            splines, scales = fit_2d_spline_profile(
                data_slice,
                alpha_slice,
                lrange=lrange,
                col_index=col_index,
                require_ngood=require_ngood,
                splinebkpt=splinebkpt,
                spaceratio=spaceratio,
                fit_scale=runsum,
            )
        else:
            splines = {}
            scales = {}
        spline_models[slnum] = splines
        spline_scales[slnum] = scales

    return spline_models, spline_scales


def just_oversample(
    flux_orig,
    alpha_orig,
    x_os,
    y_os,
    alpha_os,
    slmap,
    spline_models,
    spline_scales,
    slstart=0,
    slstop=30,
    require_ngood=8,
    trimends=False,
    pad=3,
    splinebkpt=36,
    slopelim=0.1,
    psfoptimal=False,
    oversample_factor=3,
):
    ysize, xsize = flux_orig.shape
    _, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))

    # Oversampled flux array (linear and bspline to compare)
    os_shape = (int(np.ceil(ysize * oversample_factor)), xsize)
    flux_os_linear = np.full(os_shape, np.nan)  # Linear interpolation
    flux_os_bspline_full = np.full(os_shape, np.nan)  # All bspline models
    flux_os_bspline_use = np.full(os_shape, np.nan)  # Actual bspline array applied
    flux_os_residual = np.full(os_shape, np.nan)  # Residual corrections

    # Arrays to reset with NaNs for each slice
    data_slice = np.full_like(flux_orig, np.nan)
    alpha_slice = np.full_like(flux_orig, np.nan)
    basey_slice = np.full_like(flux_orig, np.nan)

    alpha_os_slice = np.full(os_shape, np.nan)
    reset_arrays = [data_slice, basey_slice, alpha_slice, alpha_os_slice]

    # Edge limit for trimming ends
    edge_limit = int(oversample_factor)

    for slnum in range(slstart, slstop):
        # Reset holding arrays to NaN
        for reset_array in reset_arrays:
            reset_array[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = slmap == slnum
        data_slice[indx] = flux_orig[indx]
        alpha_slice[indx] = alpha_orig[indx]
        basey_slice[indx] = basey[indx]

        # Define an array that will hold all alpha values for this slice where the slope can be high
        alpha_ptsource = np.array([])

        for ii in range(xsize):
            # Are there sufficient values in this column to do anything?
            valid_data = np.isfinite(data_slice[:, ii])
            ngood = np.sum(valid_data)
            if ngood <= require_ngood:
                continue

            # Get the relevant data for this column
            col_y = basey_slice[:, ii]
            col_alpha = alpha_slice[:, ii]
            col_flux = data_slice[:, ii]

            # newy is the resampled Y pixel indices in the expanded detector frame
            # oldy is the resampled Y pixel indices in the original detector frame
            valid_y = np.isfinite(col_y)
            newy, oldy = _reindex(
                int(col_y[valid_y].min()), int(col_y[valid_y].max()), scale=oversample_factor
            )

            # Default approach is to do linear interpolation
            valid_interp = valid_y & valid_data
            interpval = np.interp(oldy, col_y[valid_interp], col_flux[valid_interp])
            if edge_limit >= 1:
                interpval[0:edge_limit] = np.nan
                interpval[-edge_limit:] = np.nan
            flux_os_linear[newy, ii] = interpval

            # Check for a spline fit for this column
            if slnum not in spline_models or ii not in spline_models[slnum]:
                continue
            spline_model = spline_models[slnum][ii]
            spline_scale = spline_scales[slnum][ii]

            # Get valid input locations and evaluate the spline
            valid_alpha = np.isfinite(col_alpha)
            col_fit = spline_model(col_alpha[valid_alpha])
            scaled_fit = col_fit * spline_scale

            # Construct the residual between spline fit and original data
            # then oversample it to output frame by linear interpolation
            residual = (col_flux[valid_alpha] - scaled_fit).astype(np.float32)
            y_interp = col_y[valid_alpha]
            valid_interp = np.isfinite(y_interp) & np.isfinite(residual)
            interpval = np.interp(oldy, y_interp[valid_interp], residual[valid_interp])
            if edge_limit >= 1:
                interpval[0:edge_limit] = np.nan
                interpval[-edge_limit:] = np.nan
            flux_os_residual[newy, ii] = interpval

            # What was the slope of the model fit prior to scaling?
            modelslope = np.abs(np.diff(col_fit, prepend=0))
            # Ensure boundaries don't look weird
            if edge_limit >= 1:
                modelslope[0:edge_limit] = 0
                modelslope[-edge_limit:] = 0

            # Where was the model slope greater than 0.1 in the native binning?
            highslope = (np.where(modelslope > slopelim))[0]

            # Add to our list of alpha values where the slope can be high for this slice
            alpha_ptsource = np.append(alpha_ptsource, col_alpha[valid_alpha][highslope])

            # Store the alpha values over sampled points in the old frame for this slice
            alpha_os_slice[newy, ii] = alpha_os[newy, ii]

            # Evaluate the bspline at the alpha for these Y locations
            oversampled_fit = spline_model(alpha_os[newy, ii]) * spline_scale
            if trimends and edge_limit >= 1:
                oversampled_fit[0:edge_limit] = np.nan
                oversampled_fit[-edge_limit:] = np.nan

            flux_os_bspline_full[newy, ii] = oversampled_fit

        # Now that our initial loop along the slice is done, we have a spline model everywhere
        # Now look at our list of alpha values where model slopes were high to figure out
        # where traces are and we actually want to use the spline model
        # Don't bother with this if not enough recorded alpha values

        if len(alpha_ptsource) > 50:
            # What is the rough native pixel size in alpha in the columns?
            native_dalpha = np.abs(np.nanmedian(np.diff(alpha_slice, axis=0)))

            avec = np.arange(splinebkpt) * native_dalpha / 2 - (native_dalpha * splinebkpt / 4)
            hist, edges = np.histogram(
                alpha_ptsource,
                bins=splinebkpt,
                range=(-native_dalpha * splinebkpt / 4, native_dalpha * splinebkpt / 4),
                density=True,
            )
            hist = hist / np.nanmax(hist)

            # Require peaks above some threshold
            peak_indices, _ = find_peaks(hist, height=0.2)
            amask = avec[peak_indices]
            for value in amask:
                indx = np.where(
                    (alpha_os_slice > value - pad * native_dalpha)
                    & (alpha_os_slice <= value + pad * native_dalpha)
                )
                flux_os_bspline_use[indx] = flux_os_bspline_full[indx]

    # Ensure that the simple interpolation didn't fill in any NaNs that it shouldn't.
    # They should be interpolated either by pixel_replace or the bspline model,
    # not this step because the simple interpolation algorithm will not be able to
    # handle them properly.
    closest_pix = (np.round(y_os).astype(int), np.round(x_os).astype(int))
    is_nan = ~np.isfinite(flux_orig[*closest_pix])
    flux_os_linear[is_nan] = np.nan

    # Insert the bspline interpolated values into the final combined oversampled array
    indx = np.where(np.isfinite(flux_os_bspline_use))
    flux_os = flux_os_linear.copy()  # Ensure we don't accidentally write into the linear data
    flux_os[indx] = flux_os_bspline_use[indx]

    # Unless we're doing a specific psf optimal extraction, add in the residual fit
    if not psfoptimal:
        log.info("Applying complex scene corrections.")
        # DRL- conflicted about this indx array
        # Using only where flux_os_bspline_use is finite will trim the slice edges a bit
        # because the spline can extend slightly beyond the linear interpolation which can
        # be bad for sources on the edge.  But requiring the residual to also be finite
        # can result in bad performance when the residual correction was really NEEDED
        # on the edge.
        indx = np.where(np.isfinite(flux_os_bspline_use) & np.isfinite(flux_os_residual))
        flux_os[indx] += flux_os_residual[indx]

    return flux_os, flux_os_bspline_full


def fit_and_oversample(
    model, lrange=50, threshsig=10, slopelim=0.1, psfoptimal=False, oversample_factor=1
):
    # slstart/stop is the slices to work on
    # pad is the padding for the replacement window
    # threshsig and slopelim can be useful to decrease if there
    # are fainter point sources that need to be spline fit too

    # todo - adapt to non-ifu data
    detector = model.meta.instrument.detector
    if detector.startswith("NRS"):
        mode = "NIRS"
        xsize, ysize = 2048, 2048
        require_ngood = 15
        splinebkpt = 62
        spaceratio = 1.6
        pad = 2
        trimends = True  # Trimming ends of the interpolation can help with bad extrapolations
        slstart, slstop = 0, 30
        chsplit = np.nan
    elif (detector == "MIRIFUSHORT") | (detector == "MIRIFULONG"):
        mode = "MIRI"
        # Note that MIRI gets rotated internally, so these are FLIPPED from usual orientation
        xsize, ysize = 1024, 1032
        require_ngood = 8
        splinebkpt = 36
        spaceratio = 1.2
        pad = 3
        trimends = (
            False  # Trimming ends is bad for MIRI, where dithers place point sources near the ends
        )

        if detector == "MIRIFUSHORT":
            slstart, slstop = 0, 38
            chsplit = 509
        else:
            slstart, slstop = 0, 28
            chsplit = 489
    else:
        raise ValueError("Unknown detector")

    basex, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))
    if mode == "NIRS":
        beta_orig = np.full((ysize, xsize), np.nan)
        alpha_orig = np.full((ysize, xsize), np.nan)
        ifu_wcs = nrs_ifu_wcs(model)
        for slice_wcs in ifu_wcs:
            x, y = gwcs.wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
            beta, alpha, _ = slice_wcs.transform("detector", "slicer", x, y)
            idx = y.astype(int), x.astype(int)
            beta_orig[*idx] = beta
            # Flip alpha so in same direction as increasing Y
            alpha_orig[*idx] = -alpha

        # NIRSpec can sort just by ordering slices by beta
        uqbeta = np.unique(beta_orig[np.isfinite(beta_orig)])
        # NIRSpec uses fluxes in normal orientation
        flux_orig = model.data
    else:
        ifu_wcs = None
        alpha_orig, beta_orig, _ = model.meta.wcs.transform(
            "detector", "alpha_beta", np.rot90(basey, k=1), np.rot90(basex, k=-1)
        )

        # todo - consider using label mapper in wcs instead
        # MIRI needs a more difficult way of sorting slices left to right
        # respecting channel boundary.
        # Hack beta of one side of the detector to ensure no overlap with the other
        beta_orig[:, 0:503] += 10
        uqbeta = []
        btemp = beta_orig[500, :]
        btemp = btemp[np.isfinite(btemp)]
        for tt in range(0, len(btemp)):
            if btemp[tt] not in uqbeta:
                uqbeta = np.append(uqbeta, btemp[tt])
        uqbeta = np.array(uqbeta)
        log.debug(f"Found {len(uqbeta)} slices.")

        # MIRI will rotate all arrays for convenience to line up with NIRSpec convention
        flux_orig = np.rot90(model.data)
        alpha_orig = np.rot90(alpha_orig)
        beta_orig = np.rot90(beta_orig)

    # Make a slice map
    slmap = np.full((ysize, xsize), np.nan)
    for ii in range(0, len(uqbeta)):
        indx = np.where(beta_orig == uqbeta[ii])
        slmap[indx] = ii

    # Set thresholding for the bspline fitting
    # Do some statistics on the overall cal file
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=AstropyUserWarning)
        overall_mean, _, overall_rms = scs(flux_orig)

    # Need to ensure that the median pixel value isn't negative, because that causes chaos
    # Subtract off that constant
    if overall_mean < 0:
        flux_orig = flux_orig - overall_mean
        overall_mean = 0

    # Define a per-slice analysis threshold (must be brighter than some level above background)
    thresh = np.full(len(uqbeta), np.nan)

    if mode == "NIRS":
        # For NIRSpec all slices have the same threshold
        thresh[:] = overall_mean + threshsig * overall_rms

    else:
        # For MIRI we need each channel to have its own threshold, particularly for Ch3/Ch4
        # since the sky is so much brighter in Ch4
        # Ch1/4
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=AstropyUserWarning)
            ch_mean, _, ch_rms = scs(flux_orig[chsplit:, :])
        slnum_in_ch = np.unique(slmap[chsplit:, :])
        for ii in range(0, len(uqbeta)):
            if ii in slnum_in_ch:
                thresh[ii] = ch_mean + threshsig * ch_rms

        # Ch2/3
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=AstropyUserWarning)
            ch_mean, _, ch_rms = scs(flux_orig[0:chsplit, :])
        slnum_in_ch = np.unique(slmap[0:chsplit, :])
        for ii in range(0, len(uqbeta)):
            if ii in slnum_in_ch:
                thresh[ii] = ch_mean + threshsig * ch_rms

    spline_models, spline_scales = just_fit(
        flux_orig,
        alpha_orig,
        slmap,
        detector,
        thresh,
        slstart,
        slstop,
        lrange,
        require_ngood,
        splinebkpt,
        spaceratio,
    )

    # If oversampling is not needed, evaluate the spline models to create the
    # profile image, store it in the model, and return
    # todo - may want to check for psf_optimal param before returning
    if oversample_factor == 1:
        profile = profile_image(flux_orig.shape, spline_models, spline_scales, slmap, alpha_orig)
        if mode == "MIRI":
            profile = np.rot90(profile, k=-1)
        model.profile = profile
        return model

    # Oversampled flux array (linear and bspline to compare)
    os_shape = (int(np.ceil(ysize * oversample_factor)), xsize)
    x_os = np.full(os_shape, np.nan)
    y_os = np.full(os_shape, np.nan)

    # Pre-compute coordinates for the new data size
    log.info("Computing oversampled coordinates")
    newy, oldy = _reindex(0, ysize - 1, scale=oversample_factor)
    y_os[:, :] = oldy[:, None]
    x_os[:, :] = basex[oldy.astype(int), :]
    if mode == "NIRS":
        alpha_os = np.full(os_shape, np.nan)
        wave_os = np.full(os_shape, np.nan)
        for slice_wcs in ifu_wcs:
            bbox = slice_wcs.bounding_box
            x_in_bounds = (x_os >= bbox[0][0]) & (x_os <= bbox[0][1])
            y_in_bounds = (y_os >= bbox[1][0]) & (y_os <= bbox[1][1])
            _, alpha, lam = slice_wcs.transform(
                "detector",
                "slicer",
                x_os[x_in_bounds & y_in_bounds],
                y_os[x_in_bounds & y_in_bounds],
            )
            alpha_os[x_in_bounds & y_in_bounds] = -alpha
            # Store wavelength, convert to um
            wave_os[x_in_bounds & y_in_bounds] = lam * 1e6
    else:
        # Because MIRI was rotated the indexing in the non-rotated frame,
        # the input coordinates need to be adjusted slightly
        alpha_os, _, wave_os = model.meta.wcs.transform(
            "detector", "alpha_beta", ysize - y_os - 1, x_os
        )

    log.info("Oversampling the flux array from the fit profile")
    flux_os, profile = just_oversample(
        flux_orig,
        alpha_orig,
        x_os,
        y_os,
        alpha_os,
        slmap,
        spline_models,
        spline_scales,
        slstart,
        slstop,
        require_ngood,
        trimends,
        pad,
        splinebkpt,
        slopelim,
        psfoptimal,
        oversample_factor,
    )

    log.info("Oversampling error and DQ arrays")
    # todo: error may need inflation to avoid underestimate later
    errors = [model.err, model.var_rnoise, model.var_poisson, model.var_flat]
    dq = model.dq
    if mode == "MIRI":
        errors = [np.rot90(err) for err in errors]
        dq = np.rot90(dq)

    # Simple linear oversample for the error arrays
    errors_os = []
    for error_array in errors:
        error_os = linear_oversample(error_array, slmap, oversample_factor, require_ngood)
        errors_os.append(error_os)

    # Nearest pixel for the dq and regions array
    # todo - update dnu to flux_estimated if no longer nan
    closest_pix = (np.round(y_os).astype(int), np.round(x_os).astype(int))
    dq_os = dq[*closest_pix]
    regions_os, pathloss_point_os, pathloss_uniform_os = None, None, None
    if model.hasattr("regions"):
        regions_os = model.regions[*closest_pix]
    elif mode == "MIRI":
        det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
        regions = np.rot90(det2ab_transform.label_mapper.mapper)
        regions_os = regions[*closest_pix]
    if model.hasattr("pathloss_point"):
        pathloss_point_os = model.pathloss_point[*closest_pix]
    if model.hasattr("pathloss_uniform"):
        pathloss_uniform_os = model.pathloss_uniform[*closest_pix]

    # Update the wcs for new pixel scale
    scale = oversample_factor
    if mode == "NIRS":
        map_pixels = Identity(1) & (Scale(1 / scale) | Shift(-(scale - 1) / (scale * 2)))

        if "coordinates" in model.meta.wcs.available_frames:
            # coordinate-based WCS
            first_transform = model.meta.wcs.pipeline[0].transform
            model.meta.wcs.pipeline[0].transform = map_pixels | first_transform
            model.meta.wcs.pipeline[0].transform.name = first_transform.name
            model.meta.wcs.pipeline[0].transform.inputs = first_transform.inputs
            model.meta.wcs.pipeline[0].transform.outputs = first_transform.outputs

            # update bounding box limits
            det2slicer_selector = model.meta.wcs.pipeline[1].transform.selector
            for slnum in range(30):
                bb = det2slicer_selector[slnum + 1].bounding_box
                bb[0], bb[1] = map_pixels.inverse(bb[0], bb[1])

        else:
            # slice-based WCS
            map_pixels &= Identity(1)
            map_pixels.name = "coord2det"
            map_pixels.inputs = ("x", "y", "name")
            map_pixels.outputs = ("x", "y", "name")
            bbox = model.meta.wcs.bounding_box
            model.meta.wcs = gwcs.WCS([("coordinates", map_pixels), *model.meta.wcs.pipeline])

            # update bounding box limits
            for slnum in range(30):
                bb = bbox[slnum]
                bb[0], bb[1], _ = map_pixels.inverse(bb[0], bb[1], slnum)
            model.meta.wcs.bounding_box = bbox
    else:
        map_pixels = (Scale(1 / scale) | Shift(-(scale - 1) / (scale * 2))) & Identity(1)
        map_pixels.name = "coord2det"
        map_pixels.inputs = ("x", "y")
        map_pixels.outputs = ("x", "y")
        model.meta.wcs = gwcs.WCS([("coordinates", map_pixels), *model.meta.wcs.pipeline])

    # If MIRI, undo all of our rotations before passing back the arrays
    if mode == "MIRI":
        flux_os = np.rot90(flux_os, k=-1)
        errors_os = [np.rot90(err, k=-1) for err in errors_os]
        dq_os = np.rot90(dq_os, k=-1)
        wave_os = np.rot90(wave_os, k=-1)
        regions_os = np.rot90(regions_os, k=-1)
        profile = np.rot90(profile, k=-1)

    model.data = flux_os
    model.err, model.var_rnoise, model.var_poisson, model.var_flat = errors_os
    model.dq = dq_os
    model.wavelength = wave_os
    model.profile = profile
    if regions_os is not None:
        model.regions = regions_os
    if pathloss_point_os is not None:
        model.pathloss_point = pathloss_point_os
    if pathloss_uniform_os is not None:
        model.pathloss_uniform = pathloss_uniform_os

    return model
