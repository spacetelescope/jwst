import logging
import warnings

import gwcs
import numpy as np
from astropy.modeling.models import Identity, Scale, Shift
from astropy.stats import sigma_clipped_stats as scs
from astropy.utils.exceptions import AstropyUserWarning
from scipy.interpolate import make_lsq_spline
from scipy.signal import find_peaks
from stdatamodels.jwst.datamodels import dqflags

from jwst.assign_wcs.nirspec import nrs_ifu_wcs
from jwst.lib.pipe_utils import match_nans_and_flags

__all__ = [
    "bspline_fit",
    "fit_2d_spline_profile",
    "linear_oversample",
    "fit_all_regions",
    "oversample_flux",
    "fit_and_oversample_ifu",
]

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
        If True, a debug log message is generated for each fitting iteration.

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
    fit_scale=None,
    lrange=50,
    col_index=None,
    require_ngood=8,
    splinebkpt=36,
    spaceratio=1.2,
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
    fit_scale : ndarray, optional
        Array of scale values to apply to the input flux before fitting.
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


def _profile_image(shape, spline_models, spline_scales, region_map, alpha):
    """
    Evaluate spline models at all pixels to generate a profile image.

    Parameters
    ----------
    shape : tuple of int
        Data shape for the output image.
    spline_models : dict of `~scipy.interpolate.BSpline`
        Spline models to evaluate.
    spline_scales : dict of float
        Scaling factors for spline models.
    region_map : ndarray
        2D image matching shape, indicating valid regions.
    alpha : ndarray
        Alpha coordinates for all pixels marked as valid regions.

    Returns
    -------
    profile : ndarray
        2D image containing the scaled spline data fit evaluated at
        the input alpha coordinates.  Values are NaN where no spline
        model was available.
    """
    profile = np.full(shape, np.nan, dtype=np.float32)
    alpha_slice = np.full(shape, np.nan, dtype=np.float32)
    profile_slice = np.full(shape, np.nan, dtype=np.float32)
    for slnum in spline_models:
        splines = spline_models[slnum]
        scales = spline_scales[slnum]
        alpha_slice[:] = np.nan
        profile_slice[:] = np.nan

        indx = region_map == slnum
        alpha_slice[indx] = alpha[indx]

        # loop over columns
        for i in range(shape[-1]):
            if i not in splines:
                continue

            # Evaluate the spline model for relevant data
            alpha_col = alpha_slice[:, i]
            valid_alpha = np.isfinite(alpha_col)
            spline_fit = scales[i] * splines[i](alpha_col[valid_alpha])
            profile_slice[:, i][valid_alpha] = spline_fit

        profile[indx] = profile_slice[indx]

    return profile


def _linear_interp(col_y, col_flux, y_interp, edge_limit, preserve_nan=True):
    """
    Perform a linear interpolation at one column.

    Parameters
    ----------
    col_y : ndarray
        Y values in the original data for the column.
    col_flux : ndarray
        Flux values in the original data for the column.
    y_interp : ndarray
        Y values to interpolate to.
    edge_limit : int
        If greater than zero, this many pixels at the edges of
        the interpolated values will be set to NaN.
    preserve_nan : bool, optional
        If True, NaNs in the input will be preserved in the output.

    Returns
    -------
    interpolated_flux : ndarray
        Interpolated flux array.
    """
    valid_data = np.isfinite(col_flux)
    valid_y = np.isfinite(col_y)

    valid_interp = valid_y & valid_data
    interpolated_flux = np.interp(y_interp, col_y[valid_interp], col_flux[valid_interp])
    if edge_limit >= 1:
        interpolated_flux[0:edge_limit] = np.nan
        interpolated_flux[-edge_limit:] = np.nan

    # Check for NaNs in the input: they should be preserved in the output
    if preserve_nan:
        closest_pix = np.round(y_interp).astype(int)
        is_nan = ~np.isfinite(col_flux[closest_pix])
        interpolated_flux[is_nan] = np.nan

    return interpolated_flux


def linear_oversample(data, region_map, oversample_factor, require_ngood, preserve_nan=True):
    """
    Oversample the input data with a linear interpolation.

    Linear interpolation is performed for each column in each region in
    the provided region map.

    Parameters
    ----------
    data : ndarray
        Original data to oversample.
    region_map : ndarray of int
        Map indicating valid regions. Values are >0 for pixels in valid
        regions, 0 otherwise.
    oversample_factor : float
        Scaling factor to oversample by.
    require_ngood : int
        Minimum number of pixels required in a column to perform an interpolation.
    preserve_nan : bool, optional
        If True, NaNs in the input will be preserved in the output.

    Returns
    -------
    os_data : ndarray
        The oversampled data array.
    """
    ysize, xsize = data.shape
    _, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))

    os_shape = (int(np.ceil(ysize * oversample_factor)), xsize)
    os_data = np.full(os_shape, np.nan, dtype=np.float32)
    edge_limit = int(oversample_factor)

    data_slice = np.full_like(data, np.nan)
    y_slice = np.full_like(data, np.nan)
    slice_numbers = np.unique(region_map[region_map > 0])
    for slnum in slice_numbers:
        data_slice[:] = np.nan
        y_slice[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = region_map == slnum
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

            os_data[newy, ii] = _linear_interp(
                col_y, col_flux, oldy, edge_limit, preserve_nan=preserve_nan
            )

    return os_data


def fit_all_regions(flux, alpha, region_map, signal_threshold, **fit_kwargs):
    """
    Fit a profile to all regions in the flux image.

    Parameters
    ----------
    flux : ndarray
        The flux image to fit.
    alpha : ndarray
        Alpha coordinates for all flux values.
    region_map : ndarray of int
        Map indicating valid regions. Values are >0 for pixels in valid
        regions, 0 otherwise.
    signal_threshold : dict of float
        Threshold values for each valid region in the region map. If
        the median peak value across columns in the region is below this
        threshold, a fit will not be attempted for that region.
    **fit_kwargs
        Keyword arguments to pass to the fitting routine (see `fit_2d_spline_profile`).

    Returns
    -------
    spline_models : dict
        Keys are region numbers, values are dicts containing a spline model for
        each column index in the region. If a spline model could not be fit, the
        column index number is not present.
    scales : dict
        Keys are region numbers, values are dicts containing a floating point scale
        for each spline model, by column index number. If a spline model could not
        be fit, the column index number is not present.
    """
    # Arrays to reset with NaNs for each slice
    data_slice = np.full_like(flux, np.nan)
    alpha_slice = np.full_like(flux, np.nan)

    spline_models = {}
    spline_scales = {}
    slice_numbers = np.unique(region_map[region_map > 0])
    for slnum in slice_numbers:
        log.info("Fitting slice %s", slnum)

        # Reset holding arrays to NaN
        data_slice[:] = np.nan
        alpha_slice[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = region_map == slnum
        data_slice[indx] = flux[indx]
        alpha_slice[indx] = alpha[indx]

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
        if medcmax > signal_threshold[slnum]:
            dospline = True

        if dospline:
            splines, scales = fit_2d_spline_profile(
                data_slice, alpha_slice, fit_scale=runsum, **fit_kwargs
            )
        else:
            splines = {}
            scales = {}
        spline_models[slnum] = splines
        spline_scales[slnum] = scales

    return spline_models, spline_scales


def oversample_flux(
    flux,
    alpha,
    region_map,
    spline_models,
    spline_scales,
    oversample_factor,
    alpha_os,
    require_ngood=8,
    trimends=False,
    pad=3,
    slopelim=0.1,
    psfoptimal=False,
):
    """
    Oversample a flux image from spline models fit to the data.

    Parameters
    ----------
    flux : ndarray
        The flux image to fit.
    alpha : ndarray
        Alpha coordinates for all flux values.
    region_map : ndarray of int
        Map indicating valid regions. Values are >0 for pixels in valid
        regions, 0 otherwise.
    spline_models : dict
        Keys are region numbers, values are dicts containing a spline model for
        each column index in the region. If a spline model could not be fit, the
        column index number is not present.
    spline_scales : dict
        Keys are region numbers, values are dicts containing a floating point scale
        for each spline model, by column index number. If a spline model could not
        be fit, the column index number is not present.
    oversample_factor : float
        Scaling factor to oversample by.
    alpha_os : ndarray
        Alpha coordinates for the oversampled array, used to evaluate spline models
        at every pixel.
    require_ngood : int, optional
        Minimum number of pixels required in a column to perform an interpolation.
    trimends : bool, optional
        If True, the edges of the evaluated spline fit will be set to NaN.
    pad : int, optional
        The number of pixels near peak data to include the spline fit for in
        the output array.
    slopelim : float, optional
        The slope limit in the normalized model fits above which the spline
        model is considered appropriate. Lower values will use spline fits
        for fainter sources.
    psfoptimal : bool, optional
        If True, residual corrections to the spline model are not included
        in the oversampled flux.

    Returns
    -------
    flux_os : ndarray
        The oversampled flux array, containing contributions from the evaluated
        spline models, linear interpolations, and residual corrections.
    profile : ndarray
        A spatial profile model, generated from the spline models evaluated at
        every pixel.
    """
    ysize, xsize = flux.shape
    _, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))

    # Oversampled flux array (linear and bspline to compare)
    os_shape = (int(np.ceil(ysize * oversample_factor)), xsize)
    flux_os_linear = np.full(os_shape, np.nan)  # Linear interpolation
    flux_os_bspline_full = np.full(os_shape, np.nan)  # All bspline models
    flux_os_bspline_use = np.full(os_shape, np.nan)  # Actual bspline array applied
    flux_os_residual = np.full(os_shape, np.nan)  # Residual corrections

    # Arrays to reset with NaNs for each slice
    data_slice = np.full_like(flux, np.nan)
    alpha_slice = np.full_like(flux, np.nan)
    basey_slice = np.full_like(flux, np.nan)

    alpha_os_slice = np.full(os_shape, np.nan)
    reset_arrays = [data_slice, basey_slice, alpha_slice, alpha_os_slice]

    # Edge limit for trimming ends
    edge_limit = int(oversample_factor)
    slice_numbers = np.unique(region_map[region_map > 0])
    splinebkpt = None
    for slnum in slice_numbers:
        # Reset holding arrays to NaN
        for reset_array in reset_arrays:
            reset_array[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = region_map == slnum
        data_slice[indx] = flux[indx]
        alpha_slice[indx] = alpha[indx]
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
            flux_os_linear[newy, ii] = _linear_interp(col_y, col_flux, oldy, edge_limit)

            # Check for a spline fit for this column
            if slnum not in spline_models or ii not in spline_models[slnum]:
                continue
            spline_model = spline_models[slnum][ii]
            spline_scale = spline_scales[slnum][ii]

            # Get the number of spline breakpoints used from the first real model
            if splinebkpt is None:
                splinebkpt = len(np.unique(spline_model.t)) - 1

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

            # Add to our list of alpha values where the slope can be high for this slice
            # and store the oversampled alpha values to check against later
            highslope = (np.where(modelslope > slopelim))[0]
            alpha_ptsource = np.append(alpha_ptsource, col_alpha[valid_alpha][highslope])
            alpha_os_slice[newy, ii] = alpha_os[newy, ii]

            # Evaluate the bspline at the oversampled alpha for this column
            oversampled_fit = spline_model(alpha_os[newy, ii]) * spline_scale
            if trimends and edge_limit >= 1:
                oversampled_fit[0:edge_limit] = np.nan
                oversampled_fit[-edge_limit:] = np.nan

            flux_os_bspline_full[newy, ii] = oversampled_fit

        # Now that our initial loop along the slice is done, we have a spline model everywhere
        # Now look at our list of alpha values where model slopes were high to figure out
        # where traces are and we actually want to use the spline model
        # Don't bother with this if not enough recorded alpha values
        if len(alpha_ptsource) > 50 and splinebkpt is not None:
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

    # Insert the bspline interpolated values into the final combined oversampled array,
    # starting from the linearly interpolated array
    flux_os = flux_os_linear
    indx = np.where(np.isfinite(flux_os_bspline_use))
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


def fit_and_oversample_ifu(
    model, threshsig=10.0, slopelim=0.1, psfoptimal=False, oversample_factor=1.0
):
    """
    Fit a spatial profile and optionally oversample an IFU datamodel.

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The input datamodel, updated in place.
    threshsig : float
        The signal threshold sigma for attempting spline fits within a slice region.
        Higher values will create spline profiles for more slices.
    slopelim : float
        The normalized slope threshold for using the spline model in oversampled
        data.  Lower values will use the spline model for fainter sources.
    psfoptimal : bool
        If True, residual corrections to the spline model are not included
        in the oversampled flux.  This option is generally appropriate for simple
        isolated point sources only.
    oversample_factor : float
        If not 1.0, then the data will be oversampled by this factor.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The datamodel, updated with a profile image and optionally oversampled
        arrays.
    """
    # todo - adapt to non-ifu data
    detector = model.meta.instrument.detector
    if detector.startswith("NRS"):
        mode = "NIRS"
        xsize, ysize = 2048, 2048
        require_ngood = 15
        splinebkpt = 62
        spaceratio = 1.6
        pad = 2
        lrange = 50

        # Trimming ends of the interpolation can help with bad extrapolations
        trimends = True
    elif (detector == "MIRIFUSHORT") | (detector == "MIRIFULONG"):
        mode = "MIRI"
        # Note that MIRI gets rotated internally, so these are FLIPPED from usual orientation
        xsize, ysize = 1024, 1032
        require_ngood = 8
        splinebkpt = 36
        spaceratio = 1.2
        pad = 3
        lrange = 50

        # Trimming ends is bad for MIRI, where dithers place point sources near the ends
        trimends = False
    else:
        raise ValueError("Unknown detector")

    basex, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))
    if mode == "NIRS":
        alpha_orig = np.full((ysize, xsize), np.nan)
        ifu_wcs = nrs_ifu_wcs(model)
        for slice_wcs in ifu_wcs:
            x, y = gwcs.wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
            _, alpha, _ = slice_wcs.transform("detector", "slicer", x, y)
            idx = y.astype(int), x.astype(int)

            # Flip alpha so in same direction as increasing Y
            alpha_orig[*idx] = -alpha

        # NIRSpec uses fluxes in normal orientation and
        # has the region map already stored in the datamodel
        flux_orig = model.data
        region_map = model.regions

        # Set up the column fitting order by detector
        if detector == "NRS1":
            # For NRS1, start on the left of detector since the tilt wrt pixels is greatest here
            col_index = range(0, xsize, 1)
        else:
            # For NRS2, start on the right of detector since the tilt wrt pixels is greatest here
            col_index = range(xsize - 1, 0, -1)
    else:
        ifu_wcs = None
        det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
        alpha_orig, _, _ = det2ab_transform(np.rot90(basey, k=1), np.rot90(basex, k=-1))

        # Region map is stored in the transform
        region_map = det2ab_transform.label_mapper.mapper

        # MIRI will rotate all arrays for convenience to line up with NIRSpec convention
        flux_orig = np.rot90(model.data)
        alpha_orig = np.rot90(alpha_orig)
        region_map = np.rot90(region_map)

        # For MIRI fitting order,  we need to start on the left and run to the middle,
        # and then on the right to the middle in order to have the middle
        # section not go too far beyond last good fit
        col_index = np.concatenate(
            [np.arange(0, xsize // 2 + 1), np.arange(xsize - 1, xsize // 2, -1)]
        )

    # Set thresholding for the bspline fitting
    # Do some statistics on the overall cal file
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=AstropyUserWarning)
        overall_mean, _, overall_rms = scs(flux_orig[region_map > 0])

    # Need to ensure that the median pixel value isn't negative, because that causes chaos
    # Subtract off that constant
    if overall_mean < 0:
        flux_orig = flux_orig - overall_mean
        overall_mean = 0

    # Define a per-slice analysis threshold (must be brighter than some level above background)
    slice_numbers = np.unique(region_map[region_map > 0])
    if mode == "NIRS":
        # For NIRSpec all slices have the same threshold
        threshold = overall_mean + threshsig * overall_rms
        signal_threshold = dict.fromkeys(slice_numbers, threshold)
    else:
        # For MIRI we need each channel to have its own threshold, particularly for Ch3/Ch4
        # since the sky is so much brighter in Ch4
        signal_threshold = dict.fromkeys(slice_numbers, np.nan)
        for channel in [100, 200, 300, 400]:
            ch_data = (region_map >= channel) & (region_map < channel + 100)
            if not np.any(ch_data):
                continue
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=AstropyUserWarning)
                ch_mean, _, ch_rms = scs(flux_orig[ch_data])
            for slnum in slice_numbers:
                if channel <= slnum < channel + 100:
                    signal_threshold[slnum] = ch_mean + threshsig * ch_rms

    fit_kwargs = {
        "lrange": lrange,
        "col_index": col_index,
        "require_ngood": require_ngood,
        "splinebkpt": splinebkpt,
        "spaceratio": spaceratio,
    }
    spline_models, spline_scales = fit_all_regions(
        flux_orig, alpha_orig, region_map, signal_threshold, **fit_kwargs
    )

    # If oversampling is not needed, evaluate the spline models to create the
    # profile image, store it in the model, and return
    # todo - may want to check for psf_optimal param before returning
    if oversample_factor == 1:
        profile = _profile_image(
            flux_orig.shape, spline_models, spline_scales, region_map, alpha_orig
        )
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
    flux_os, profile = oversample_flux(
        flux_orig,
        alpha_orig,
        region_map,
        spline_models,
        spline_scales,
        oversample_factor,
        alpha_os,
        require_ngood=require_ngood,
        trimends=trimends,
        pad=pad,
        slopelim=slopelim,
        psfoptimal=psfoptimal,
    )

    log.info("Oversampling error and DQ arrays")
    # todo: error may need inflation to avoid underestimate later
    errors = [model.err, model.var_rnoise, model.var_poisson, model.var_flat]
    dq = model.dq
    if mode == "MIRI":
        errors = [np.rot90(err) for err in errors]
        dq = np.rot90(dq)

    # Nearest pixel interpolation for the dq and regions array
    closest_pix = (np.round(y_os).astype(int), np.round(x_os).astype(int))
    dq_os = dq[*closest_pix]
    regions_os = region_map[*closest_pix]
    pathloss_point_os, pathloss_uniform_os = None, None
    if model.hasattr("pathloss_point"):
        pathloss_point_os = model.pathloss_point[*closest_pix]
    if model.hasattr("pathloss_uniform"):
        pathloss_uniform_os = model.pathloss_uniform[*closest_pix]

    # Update the DQ image for pixels that used to be NaN, now replaced by spline interpolation.
    # Remove the DO_NOT_USE flag, add FLUX_ESTIMATED
    is_estimated = ~np.isnan(flux_os) & ((dq_os & dqflags.pixel["DO_NOT_USE"]) > 0)
    dq_os[is_estimated] ^= dqflags.pixel["DO_NOT_USE"]
    dq_os[is_estimated] |= dqflags.pixel["FLUX_ESTIMATED"]

    # Simple linear oversample for the error arrays
    errors_os = []
    for error_array in errors:
        error_os = linear_oversample(
            error_array, region_map, oversample_factor, require_ngood, preserve_nan=False
        )

        # Restore NaNs from the input, except at the estimated locations
        is_nan = ~np.isfinite(error_array[closest_pix])
        error_os[is_nan & ~is_estimated] = np.nan

        errors_os.append(error_os)

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

    # Make sure NaNs and DO_NOT_USE flags match in all extensions
    match_nans_and_flags(model)

    return model
