import logging
import warnings

import gwcs
import numpy as np
from astropy.modeling.models import Identity, Scale, Shift
from astropy.stats import sigma_clipped_stats as scs
from astropy.utils.exceptions import AstropyUserWarning
from scipy.signal import find_peaks
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.adaptive_trace_model.bspline import bspline_fit
from jwst.assign_wcs.nirspec import nrs_ifu_wcs
from jwst.lib.pipe_utils import match_nans_and_flags

__all__ = [
    "fit_2d_spline_trace",
    "linear_oversample",
    "fit_all_regions",
    "oversample_flux",
    "fit_and_oversample",
]

log = logging.getLogger(__name__)


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


def fit_2d_spline_trace(
    flux,
    alpha,
    fit_scale=None,
    lrange=50,
    col_index=None,
    require_ngood=10,
    spline_bkpt=50,
    space_ratio=1.2,
):
    """
    Create a trace model from spline fits to a single slit/slice image.

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
    spline_bkpt : int, optional
        Number of spline breakpoints (knots).
    space_ratio : float, optional
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
        col_alpha = alpha[:, i]
        ngood = np.sum(np.isfinite(col_flux) & np.isfinite(col_alpha))
        if ngood <= require_ngood:
            continue

        # Get local alpha and flux values for fitting
        lstart = np.max([i - lrange, 0])
        lstop = np.min([i + lrange, xsize])
        local_alpha = alpha[:, lstart:lstop]
        local_data = scaled_flux[:, lstart:lstop]

        # Trim to finite values
        finite_values = np.isfinite(local_alpha) & np.isfinite(local_data)
        local_alpha = local_alpha[finite_values]
        local_data = local_data[finite_values]

        # Sort by alpha
        idx = np.argsort(local_alpha)
        local_alpha = local_alpha[idx]
        local_data = local_data[idx]

        # Fit a bspline to the local data
        try:
            bspline = bspline_fit(
                local_alpha,
                local_data,
                nbkpts=spline_bkpt,
                wrapsig_low=2.5,
                wrapsig_high=2.5,
                wrapiter=3,
                space_ratio=space_ratio,
                verbose=False,
            )

            # If this routine could not get a fit (returned None) use the saved fit
            if bspline is None:
                spline_model = spline_model_save
            else:
                spline_model = bspline
        except (ValueError, RuntimeError) as err:
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
        idx = np.isfinite(col_alpha)
        col_alpha = col_alpha[idx]
        col_flux = col_flux[idx]
        col_fit = spline_model(col_alpha)

        # Determine the normalization by the weighted mean ratio between model and data
        # Weights are based on the model so that we can reject outliers
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)

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

    For example, with oversample scale = 2, [0, 1, 2] goes to
    old_x = [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25], for
    new_x = [0, 1, 2, 3, 4, 5].
    With oversample scale = 3, [0, 1, 2] goes to
    old_x = [-0.33, 0, 0.33, 0.67, 1, 1.33, 1.67, 2, 2.33], for
    new_x = [0, 1, 2, 3, 4, 5, 6, 7, 8].

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


def _is_compact_source(
    alpha_slice, alpha_ptsource, native_dalpha, spline_bkpt, pad=3, require_npt=50
):
    """
    Determine which pixels within a slice contain a compact source.

    Parameters
    ----------
    alpha_slice : ndarray
        Alpha coordinates for the output slice. If oversampling is performed,
        these should be the oversampled coordinates.
    alpha_ptsource : ndarray
        Array of alpha values for modeled flux that met the slope limit threshold.
    native_dalpha : float
        Approximate native pixel size in alpha, along the columns.
    spline_bkpt : int
        The number of breakpoints used in the spline modeling.
    pad : int, optional
        The number of pixels near peak data to include the spline fit for in
        the output array.
    require_npt : int, optional
        The minimum required number of high-slope data points to consider any
        pixels to be compact.

    Returns
    -------
    is_compact : ndarray
        Boolean array matching the shape of ``alpha_slice``, where True
        indicates a pixel containing a compact source.
    """
    is_compact = np.full(alpha_slice.shape, False)

    # If there is not enough high slope data or no spline models were found,
    # just return False for all data
    if len(alpha_ptsource) < require_npt or spline_bkpt is None:
        return is_compact

    # Bin the alpha coordinates for the high slope locations
    avec = np.arange(spline_bkpt) * native_dalpha / 2 - (native_dalpha * spline_bkpt / 4)
    hist, edges = np.histogram(
        alpha_ptsource,
        bins=spline_bkpt,
        range=(-native_dalpha * spline_bkpt / 4, native_dalpha * spline_bkpt / 4),
        density=True,
    )
    hist = hist / np.nanmax(hist)

    # Require peaks above some threshold
    peak_indices, _ = find_peaks(hist, height=0.2)
    amask = avec[peak_indices]

    # Flag regions near the compact source, with some padding
    for value in amask:
        indx = (alpha_slice > value - pad * native_dalpha) & (
            alpha_slice <= value + pad * native_dalpha
        )
        is_compact[indx] = True
    return is_compact


def _trace_image(shape, spline_models, spline_scales, region_map, alpha, slope_limit=0.1, pad=3):
    """
    Evaluate spline models at all pixels to generate a trace image.

    The trace image will be NaN wherever a spline model was not fit and
    wherever the source is not compact enough for the spline model
    to be appropriate. The ``slope_limit`` parameter controls the decision
    for compact source regions.

    Parameters
    ----------
    shape : tuple of int
        Data shape for the output image.
    spline_models : dict
        Spline models to evaluate.
    spline_scales : dict
        Scaling factors for spline models.
    region_map : ndarray
        2D image matching shape, mapping valid region numbers.
    alpha : ndarray
        Alpha coordinates for all pixels marked as valid regions.
    slope_limit : float, optional
        The slope limit in the normalized model fits above which the spline
        model is considered appropriate. Lower values will use spline fits
        for fainter sources. If less than or equal to zero, the spline fits
        will always be used.
    pad : int, optional
        The number of pixels near peak data to include the spline fit for in
        the output array.

    Returns
    -------
    trace_used : ndarray
        2D image containing the scaled spline data fit evaluated at
        the input alpha coordinates for compact source regions.  Values are
        NaN where no spline model was available and where the source
        was below the slope limit.
    full_trace: ndarray
        2D image containing the scaled spline data fit evaluated at
        all pixels. Values are NaN where no spline model was available.
    """
    trace_used = np.full(shape, np.nan, dtype=np.float32)
    full_trace = np.full(shape, np.nan, dtype=np.float32)
    alpha_slice = np.full(shape, np.nan, dtype=np.float32)
    trace_slice = np.full(shape, np.nan, dtype=np.float32)
    spline_bkpt = None
    for slnum in spline_models:
        splines = spline_models[slnum]
        scales = spline_scales[slnum]
        alpha_slice[:] = np.nan
        trace_slice[:] = np.nan

        indx = region_map == slnum
        alpha_slice[indx] = alpha[indx]

        # Define a list that will hold all alpha values for this
        # slice where the slope is high
        alpha_ptsource = []

        # loop over columns
        for i in range(shape[-1]):
            if i not in splines:
                continue

            # Evaluate the spline model for relevant data
            col_alpha = alpha_slice[:, i]
            valid_alpha = np.isfinite(col_alpha)
            col_fit = splines[i](col_alpha[valid_alpha])

            # Set the edges to NaN to avoid edge effects
            col_fit[0] = np.nan
            col_fit[-1] = np.nan

            scaled_fit = scales[i] * col_fit
            trace_slice[:, i][valid_alpha] = scaled_fit

            # Get the slope of the model fit prior to scaling
            model_slope = np.abs(np.diff(col_fit, prepend=0))

            # Ensure boundaries don't look weird
            model_slope[0] = 0
            model_slope[-1] = 0

            highslope = (np.where(model_slope > slope_limit))[0]
            alpha_ptsource.append(col_alpha[valid_alpha][highslope])

            # Get the number of spline breakpoints used from the first real model
            if spline_bkpt is None:
                spline_bkpt = len(np.unique(splines[i].t)) - 1

        full_trace[indx] = trace_slice[indx]
        if slope_limit <= 0:
            # Always use the spline fit in this case
            trace_used[indx] = trace_slice[indx]
        else:
            if len(alpha_ptsource) > 0:
                alpha_ptsource = np.concatenate(alpha_ptsource)

            native_dalpha = np.abs(np.nanmedian(np.diff(alpha_slice, axis=0)))
            compact_locations = _is_compact_source(
                alpha_slice, alpha_ptsource, native_dalpha, spline_bkpt, pad
            )
            trace_used[compact_locations] = trace_slice[compact_locations]

            total_used = np.sum(compact_locations)
            log.debug(
                f"Using {total_used}/{np.sum(indx)} pixels from the spline model for slice {slnum}"
            )

    return trace_used, full_trace


def _linear_interp(col_y, col_flux, y_interp, edge_limit=0, preserve_nan=True):
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
    edge_limit : int, optional
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


def linear_oversample(
    data, region_map, oversample_factor, require_ngood, edge_limit=0, preserve_nan=True
):
    """
    Oversample the input data with a linear interpolation.

    Linear interpolation is performed for each column in each region in
    the provided region map.

    Parameters
    ----------
    data : ndarray
        Original data to oversample.
    region_map : ndarray of int
        Map containing the slice or slit number for valid regions.
        Values are >0 for pixels in valid regions, 0 otherwise.
    oversample_factor : float
        Scaling factor to oversample by.
    require_ngood : int
        Minimum number of pixels required in a column to perform an interpolation.
    edge_limit : int, optional
        If greater than zero, this many pixels at the edges of
        the interpolated values will be set to NaN.
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
                col_y, col_flux, oldy, edge_limit=edge_limit, preserve_nan=preserve_nan
            )

    return os_data


def fit_all_regions(flux, alpha, region_map, signal_threshold, **fit_kwargs):
    """
    Fit a trace model to all regions in the flux image.

    Parameters
    ----------
    flux : ndarray
        The flux image to fit.
    alpha : ndarray
        Alpha coordinates for all flux values.
    region_map : ndarray of int
        Map containing the slice or slit number for valid regions.
        Values are >0 for pixels in valid regions, 0 otherwise.
    signal_threshold : dict
        Threshold values for each valid region in the region map. If
        the median peak value across columns in the region is below this
        threshold, a fit will not be attempted for that region.
    **fit_kwargs
        Keyword arguments to pass to the fitting routine (see `fit_2d_spline_trace`).

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
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            collapse = np.nanmax(data_slice, axis=0)

        # Median column max across all columns
        medcmax = np.nanmedian(collapse)

        # Is medcmax over threshold?  If so, do bspline for this slice.
        dospline = False
        if medcmax > signal_threshold[slnum]:
            dospline = True

        if dospline:
            splines, scales = fit_2d_spline_trace(
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
    require_ngood=10,
    slope_limit=0.1,
    psf_optimal=False,
    trim_ends=False,
    pad=3,
):
    """
    Oversample a flux image from spline models fit to the data.

    For each column in each slice or slit in the region map:

    1. Check if there are enough valid data points to proceed.
    2. Compute oversampled coordinates corresponding to the input column.
    3. Linearly interpolate flux values onto the oversampled column.
    4. If a spline fit is available, evaluate it for the original column
       coordinates.
    5. Construct a residual between the spline fit and the original column.
       data, then linearly interpolate the residual onto the oversampled
       column.
    6. Compute the slope of each column pixel as the difference between the
       normalized spline model at that pixel and its immediate neighbor.
    7. Evaluate the spline model at the oversampled coordinates.

    The oversampled flux for each slice or slit is set from the spline flux
    plus the interpolated residual, for pixels where the slope exceeds the
    ``slope_limit``.  Otherwise, the flux is set to the linearly interpolated
    value.

    Parameters
    ----------
    flux : ndarray
        The flux image to fit.
    alpha : ndarray
        Alpha coordinates for all flux values.
    region_map : ndarray of int
        Map containing the slice or slit number for valid regions.
        Values are >0 for pixels in valid regions, 0 otherwise.
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
    slope_limit : float, optional
        The slope limit in the normalized model fits above which the spline
        model is considered appropriate. Lower values will use spline fits
        for fainter sources. If less than or equal to zero, the spline fits
        will always be used.
    psf_optimal : bool, optional
        If True, residual corrections to the spline model are not included
        in the oversampled flux.
    trim_ends : bool, optional
        If True, the edges of the evaluated spline fit will be set to NaN.
    pad : int, optional
        The number of pixels near peak data to include the spline fit for in
        the output array.

    Returns
    -------
    flux_os : ndarray
        The oversampled flux array, containing contributions from the evaluated
        spline models, linear interpolations, and residual corrections.
    trace_used : ndarray
        A trace model, generated from the spline models evaluated at
        pixels containing a compact source.
    full_trace : ndarray
        A trace model, generated from the spline models evaluated at
        every pixel.
    linear_flux : ndarray
        The flux linearly interpolated onto the oversampled grid.
    residual_flux : ndarray
        Residuals between the spline modeled data and the original flux,
        linearly interpolated onto the oversampled grid.
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
    spline_bkpt = None
    for slnum in slice_numbers:
        # Reset holding arrays to NaN
        for reset_array in reset_arrays:
            reset_array[:] = np.nan

        # Copy the relevant data for this slice into the holding arrays
        indx = region_map == slnum
        data_slice[indx] = flux[indx]
        alpha_slice[indx] = alpha[indx]
        basey_slice[indx] = basey[indx]

        # Define a list that will hold all alpha values for this slice
        # where the slope is high
        alpha_ptsource = []

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
            flux_os_linear[newy, ii] = _linear_interp(col_y, col_flux, oldy, edge_limit=edge_limit)

            # Check for a spline fit for this column
            if slnum not in spline_models or ii not in spline_models[slnum]:
                continue
            spline_model = spline_models[slnum][ii]
            spline_scale = spline_scales[slnum][ii]

            # Get the number of spline breakpoints used from the first real model
            if spline_bkpt is None:
                spline_bkpt = len(np.unique(spline_model.t)) - 1

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
            model_slope = np.abs(np.diff(col_fit, prepend=0))
            # Ensure boundaries don't look weird
            if edge_limit >= 1:
                model_slope[0:edge_limit] = 0
                model_slope[-edge_limit:] = 0

            # Add to our list of alpha values where the slope can be high for this slice
            highslope = (np.where(model_slope > slope_limit))[0]
            alpha_ptsource.append(col_alpha[valid_alpha][highslope])

            # Store the oversampled alpha values to check against later
            alpha_os_slice[newy, ii] = alpha_os[newy, ii]

            # Evaluate the bspline at the oversampled alpha for this column
            oversampled_fit = spline_model(alpha_os[newy, ii]) * spline_scale
            if trim_ends and edge_limit >= 1:
                oversampled_fit[0:edge_limit] = np.nan
                oversampled_fit[-edge_limit:] = np.nan

            flux_os_bspline_full[newy, ii] = oversampled_fit

        # Now that our initial loop along the slice is done, we have a spline model everywhere

        # Now look at our list of alpha values where model slopes were high to figure out
        # where traces are and we actually want to use the spline model
        if slope_limit <= 0:
            # Always use the spline fit in this case
            flux_os_bspline_use = flux_os_bspline_full
        else:
            if len(alpha_ptsource) > 0:
                alpha_ptsource = np.concatenate(alpha_ptsource)
            native_dalpha = np.abs(np.nanmedian(np.diff(alpha_slice, axis=0)))
            compact_locations = _is_compact_source(
                alpha_os_slice, alpha_ptsource, native_dalpha, spline_bkpt, pad
            )
            flux_os_bspline_use[compact_locations] = flux_os_bspline_full[compact_locations]

            total_used = np.sum(compact_locations)
            log.debug(
                f"Using {total_used}/{np.sum(indx)} pixels from the spline model for slice {slnum}"
            )

    # Insert the bspline interpolated values into the final combined oversampled array,
    # starting from the linearly interpolated array
    flux_os = flux_os_linear
    indx = np.where(np.isfinite(flux_os_bspline_use))
    flux_os[indx] = flux_os_bspline_use[indx]

    # Unless we're doing a specific psf optimal extraction, add in the residual fit
    if not psf_optimal:
        log.info("Applying complex scene corrections.")
        # DRL- conflicted about this indx array
        # Using only where flux_os_bspline_use is finite will trim the slice edges a bit
        # because the spline can extend slightly beyond the linear interpolation which can
        # be bad for sources on the edge.  But requiring the residual to also be finite
        # can result in bad performance when the residual correction was really NEEDED
        # on the edge.
        indx = np.where(np.isfinite(flux_os_bspline_use) & np.isfinite(flux_os_residual))
        flux_os[indx] += flux_os_residual[indx]

    return flux_os, flux_os_bspline_use, flux_os_bspline_full, flux_os_linear, flux_os_residual


def _set_fit_kwargs(detector, xsize):
    """
    Set optional parameters for spline fits by detector.

    Parameters
    ----------
    detector : str
        Detector name.
    xsize : int
        Input size for the data, along the dispersion axis. Used
        to determine the column index order for spline fits.

    Returns
    -------
    fit_kwargs : dict
        Optional parameter settings to pass to the ``fit_all_regions``
        function.

    Raises
    ------
    ValueError
        If the input detector is not supported.
    """
    # Empirical parameters for this mode
    if detector.startswith("NRS"):
        require_ngood = 15
        spline_bkpt = 62
        lrange = 50

        # This factor of 1.6 was dialed based on inspection of the results
        # as sampling gets progressively worse for NIRSpec detectors
        space_ratio = 1.6

        # Set up the column fitting order by detector
        if detector == "NRS1":
            # For NRS1, start on the left of detector since the tilt wrt pixels is greatest here
            col_index = range(0, xsize, 1)
        else:
            # For NRS2, start on the right of detector since the tilt wrt pixels is greatest here
            col_index = range(xsize - 1, -1, -1)

    elif detector.startswith("MIR"):
        require_ngood = 8
        spline_bkpt = 36
        lrange = 50
        space_ratio = 1.2

        # For MIRI fitting order,  we need to start on the left and run to the middle,
        # and then on the right to the middle in order to have the middle
        # section not go too far beyond last good fit
        col_index = np.concatenate(
            [np.arange(0, xsize // 2 + 1), np.arange(xsize - 1, xsize // 2, -1)]
        )
    else:
        raise ValueError("Unknown detector")

    fit_kwargs = {
        "lrange": lrange,
        "col_index": col_index,
        "require_ngood": require_ngood,
        "spline_bkpt": spline_bkpt,
        "space_ratio": space_ratio,
    }
    return fit_kwargs


def _set_oversample_kwargs(detector):
    """
    Set optional parameters for oversampling by detector.

    Parameters
    ----------
    detector : str
        Detector name.

    Returns
    -------
    oversample_kwargs : dict
        Optional parameter settings to pass to the ``oversample_flux``
        function.

    Raises
    ------
    ValueError
        If the input detector is not supported.
    """
    if detector.startswith("NRS"):
        # Trimming ends of the interpolation can help with bad extrapolations
        pad = 2
        trim_ends = True
    elif detector.startswith("MIR"):
        # Trimming ends is bad for MIRI, where dithers place point sources near the ends
        pad = 3
        trim_ends = False
    else:
        raise ValueError("Unknown detector")

    oversample_kwargs = {"pad": pad, "trim_ends": trim_ends}
    return oversample_kwargs


def _get_alpha_nrs_ifu(ifu_wcs, xsize, ysize):
    """
    Get alpha coordinates for NIRSpec IFU corresponding to the original data array.

    Parameters
    ----------
    ifu_wcs : list of `~gwcs.WCS`
        List of WCS objects, one per slice.
    xsize : int
        X-size for the data array.
    ysize : int
        Y-size for the data array.

    Returns
    -------
    alpha_orig : ndarray
        Alpha coordinates for the data array, with shape (ysize, xsize).
    """
    alpha_orig = np.full((ysize, xsize), np.nan)
    for slice_wcs in ifu_wcs:
        x, y = gwcs.wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
        _, alpha, _ = slice_wcs.transform("detector", "slicer", x, y)
        idx = y.astype(int), x.astype(int)

        # Flip alpha so in same direction as increasing Y
        alpha_orig[*idx] = -alpha
    return alpha_orig


def _get_alpha_mir_mrs(wcs, xsize, ysize):
    """
    Get alpha coordinates for MIRI MRS corresponding to the original data array.

    Parameters
    ----------
    wcs : `~gwcs.WCS`
        WCS object.
    xsize : int
        X-size for the data array.
    ysize : int
        Y-size for the data array.

    Returns
    -------
    alpha_orig : ndarray
        Alpha coordinates for the data array, with shape (ysize, xsize).
    """
    x, y = np.meshgrid(np.arange(xsize), np.arange(ysize))
    det2ab = wcs.get_transform("detector", "alpha_beta")
    alpha_orig, _, _ = det2ab(x, y)
    return alpha_orig


def _get_oversampled_coords_nrs_ifu(ifu_wcs, x_os, y_os):
    """
    Get alpha coordinates for NIRSpec IFU corresponding to the oversampled data array.

    Parameters
    ----------
    ifu_wcs : list of `~gwcs.WCS`
        List of WCS objects, one per slice.
    x_os : int
        X-size for the oversampled data array.
    y_os : int
        Y-size for the oversampled data array.

    Returns
    -------
    alpha_os : ndarray
        Alpha coordinates for the data array, with shape (y_os, x_os).
    wave_os : ndarray
        Wavelength coordinates for the data array, with shape (y_os, x_os),
        in um.
    """
    os_shape = x_os.shape
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

    return alpha_os, wave_os


def _inflate_error(error_array, extname, oversample_factor):
    """
    Inflate error or variance arrays to account for oversampling.

    Errors are increased by a factor dependent on the oversampling ratio
    in order to account for the covariance introduced by the oversampling.

    The inflation factor was determined empirically for IFU data
    by comparing the reported error of single-spaxel spectra and aperture-summed
    spectra, following ``cube_build`` on an oversampled image.

    Empirically, based on the RMS of the aperture-summed spectrum in a line-free
    region of  a stellar spectrum, the true SNR does not change much (< 4%) between
    N=1 and N=2/3/4. In contrast the reported SNR increases by an amount well fit
    by X = 0.23N + 0.77. I.e., X=1 for N=1, and X=1.46 for N=3.  This does not account
    for variations in individual pixels, but to first order, inflating by this X
    factor when the oversampling is performed will produce data cubes in which the
    SNR is mostly preserved accurately.  Per-pixel errors in the oversampled
    product are not accurately reported by the inflated errors, but the oversampled
    product should be considered primarily an intermediate data product; the
    errors in the resampled cube are more important.

    Parameters
    ----------
    error_array : ndarray
        Error or variance image to inflate. Updated in place.
    extname : {"err", "var_rnoise", "var_poisson", "var_flat"}
        Extension name.
    oversample_factor : float
        The oversampling factor used.
    """
    inflation_factor = 0.23 * oversample_factor + 0.77
    if str(extname).lower().startswith("var"):
        error_array *= inflation_factor**2
    else:
        error_array *= inflation_factor


def _update_wcs_nrs_ifu(wcs, map_pixels):
    """
    Update a NIRSpec IFU WCS to include the oversampling transform.

    Parameters
    ----------
    wcs : `~gwcs.WCS`
        The WCS object, including transforms for all slices.
        May be either coordinate-based or slice-based.
    map_pixels : `~astropy.modeling.models.Model`
        Model that transforms from oversampled pixels to original detector
        pixels, to be prepended to the WCS pipeline.

    Returns
    -------
    wcs : `~gwcs.WCS`
        The updated WCS.  If the input WCS was coordinate-based,
        then the new transform is prepended to the existing "coordinates"
        transform.  If it was slice-based, a new WCS pipeline is created
        with "coordinates" as the input frame, containing the new transform.
    """
    if "coordinates" in wcs.available_frames:
        # coordinate-based WCS: update the existing transform with the new mapping
        first_transform = wcs.pipeline[0].transform
        wcs.pipeline[0].transform = map_pixels | first_transform
        wcs.pipeline[0].transform.name = first_transform.name
        wcs.pipeline[0].transform.inputs = first_transform.inputs
        wcs.pipeline[0].transform.outputs = first_transform.outputs

        # update bounding box limits in place
        det2slicer_selector = wcs.pipeline[1].transform.selector
        for slnum in range(30):
            bb = det2slicer_selector[slnum + 1].bounding_box
            bb[0], bb[1] = map_pixels.inverse(bb[0], bb[1])

    else:
        # slice-based WCS
        map_pixels &= Identity(1)
        map_pixels.name = "coord2det"
        map_pixels.inputs = ("x", "y", "name")
        map_pixels.outputs = ("x", "y", "name")
        bbox = wcs.bounding_box
        frame = gwcs.coordinate_frames.Frame2D(name="coordinates", axes_order=(0, 1))
        wcs = gwcs.WCS([(frame, map_pixels), *wcs.pipeline])

        # update bounding box limits
        for slnum in range(30):
            bb = bbox[slnum]
            bb[0], bb[1], _ = map_pixels.inverse(bb[0], bb[1], slnum)
        wcs.bounding_box = bbox
    return wcs


def _update_wcs(wcs, map_pixels):
    """
    Update a WCS to include the oversampling transform.

    Appropriate to the MIRI MRS WCS or slit-like WCS objects, following ``extract_2d``.

    Parameters
    ----------
    wcs : `~gwcs.WCS`
        The WCS object, including transforms for all slices.
    map_pixels : `~astropy.modeling.models.Model`
        Model that transforms from oversampled pixels to original detector
        pixels, to be prepended to the WCS pipeline.

    Returns
    -------
    wcs : `~gwcs.WCS`
        A new WCS pipeline, with "coordinates" as the input frame, containing the
        new transform.
    """
    map_pixels.name = "coord2det"
    map_pixels.inputs = ("x", "y")
    map_pixels.outputs = ("x", "y")
    frame = gwcs.coordinate_frames.Frame2D(name="coordinates", axes_order=(0, 1))
    new_wcs = gwcs.WCS([(frame, map_pixels), *wcs.pipeline])
    return new_wcs


def _intermediate_models(model, data_arrays):
    """
    Make new datamodels for intermediate data arrays.

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The input datamodel. Metadata will be copied from it.
    data_arrays : list of ndarray or None
        Data arrays to save.  If None, the model returned is also None.

    Returns
    -------
    new_models : list of `~stdatamodels.jwst.datamodels.IFUImageModel` or None
        A list of datamodels containing the input data arrays.
    """
    new_models = []
    for data in data_arrays:
        if data is None:
            new_model = None
        else:
            new_model = datamodels.IFUImageModel(data)
            new_model.update(model)
        new_models.append(new_model)
    return new_models


def fit_and_oversample(
    model,
    fit_threshold=10.0,
    slope_limit=0.1,
    psf_optimal=False,
    oversample_factor=1.0,
    return_intermediate_models=False,
):
    """
    Fit a trace model and optionally oversample an IFU datamodel.

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The input datamodel, updated in place.
    fit_threshold : float, optional
        The signal threshold sigma for attempting spline fits within a slice region.
        Lower values will create spline traces for more slices.  If less than or
        equal to 0, all slices will be fit.
    slope_limit : float, optional
        The normalized slope threshold for using the spline model in oversampled
        data.  Lower values will use the spline model for fainter sources. If less
        than or equal to 0, the spline model will always be used.
    psf_optimal : bool, optional
        If True, residual corrections to the spline model are not included
        in the oversampled flux.  This option is generally appropriate for simple
        isolated point sources only.  If set, ``slope_limit`` and ``fit_threshold``
        values are ignored and spline fits are attempted and used for all data.
    oversample_factor : float, optional
        If not 1.0, then the data will be oversampled by this factor.
    return_intermediate_models : bool, optional
        If True, additional image models will be returned, containing the full
        spline model, the spline model as used for compact sources, the residual
        model, and the linearly interpolated data.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The datamodel, updated with a trace image and optionally oversampled
        arrays.
    full_spline_model : `~stdatamodels.jwst.datamodels.IFUImageModel`, optional
        The spline model evaluated at all pixels. Returned only if
        ``return_intermediate_models`` is True.
    source_spline_model : `~stdatamodels.jwst.datamodels.IFUImageModel`, optional
        The spline model evaluated at compact source locations only.
        Returned only if ``return_intermediate_models`` is True.
    linear_model : `~stdatamodels.jwst.datamodels.IFUImageModel` or None, optional
        All data linearly interpolated onto the oversampled grid
        Returned only if ``return_intermediate_models`` is True.  Will be None if
        ``oversample_factor`` is 1.0.
    residual_model : `~stdatamodels.jwst.datamodels.IFUImageModel` or None, optional
        Residuals from the spline fit, linearly interpolated onto the oversampled grid
        Returned only if ``return_intermediate_models`` is True.  Will be None if
        ``oversample_factor`` is 1.0.
    """
    # Check parameters
    if psf_optimal:
        log.info("Ignoring fit threshold and slope limit for psf_optimal=True")
        fit_threshold = 0
        slope_limit = 0

    # Get input data coordinates
    detector = model.meta.instrument.detector
    ysize, xsize = model.data.shape
    if detector.startswith("NRS"):
        rotate = False
        if isinstance(model, datamodels.IFUImageModel):
            mode = "NRS_IFU"
            wcs = nrs_ifu_wcs(model)
            alpha_orig = _get_alpha_nrs_ifu(wcs, xsize, ysize)

            # the region map is already stored in the datamodel
            region_map = model.regions
        else:
            raise ValueError("Unsupported mode")

    elif detector.startswith("MIR"):
        rotate = True
        if isinstance(model, datamodels.IFUImageModel):
            mode = "MIR_MRS"
            wcs = model.meta.wcs
            alpha_orig = _get_alpha_mir_mrs(wcs, xsize, ysize)

            # Region map is stored in the transform
            det2ab_transform = wcs.get_transform("detector", "alpha_beta")
            region_map = det2ab_transform.label_mapper.mapper.copy()
        else:
            raise ValueError("Unsupported mode")
    else:
        raise ValueError("Unknown detector")

    # Rotate input data if needed
    flux_orig = model.data
    if rotate:
        xsize, ysize = ysize, xsize
        flux_orig = np.rot90(flux_orig)
        alpha_orig = np.rot90(alpha_orig)
        region_map = np.rot90(region_map)

    # Set thresholding for the bspline fitting
    # Do some statistics on the overall cal file
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=AstropyUserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        overall_mean, _, overall_rms = scs(flux_orig[region_map > 0])
    overall_mean = 0 if ~np.isfinite(overall_mean) else overall_mean
    overall_rms = 0 if ~np.isfinite(overall_rms) else overall_rms

    # Need to ensure that the median pixel value isn't negative, because that causes chaos
    # Subtract off that constant
    if overall_mean < 0:
        flux_orig = flux_orig - overall_mean
        overall_mean = 0

    # Define a per-slice analysis threshold (must be brighter than some level above background)
    slice_numbers = np.unique(region_map[region_map > 0])
    if fit_threshold <= 0:
        # In this case, all slices should be fit, so make the threshold
        # lower than any real signal
        signal_threshold = dict.fromkeys(slice_numbers, -np.inf)
    else:
        if mode == "MIR_MRS":
            # For MIRI MRS we need each channel to have its own threshold, particularly
            # for Ch3/Ch4 since the sky is so much brighter in Ch4
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
                        signal_threshold[slnum] = ch_mean + fit_threshold * ch_rms
        else:
            # For NIRSpec IFU, all regions have the same threshold
            threshold = overall_mean + fit_threshold * overall_rms
            signal_threshold = dict.fromkeys(slice_numbers, threshold)

    # Fit spline models to all regions
    fit_kwargs = _set_fit_kwargs(detector, xsize)
    spline_models, spline_scales = fit_all_regions(
        flux_orig, alpha_orig, region_map, signal_threshold, **fit_kwargs
    )

    # If oversampling is not needed, evaluate the spline models to create the
    # trace image, store it in the model, and return.
    # In the future, it might be useful to update the SCI extension here for the
    # psf_optimal=True case, even when oversample=1, but for now, we will leave
    # data unmodified.
    oversample_kwargs = _set_oversample_kwargs(detector)
    if oversample_factor == 1:
        trace_used, full_trace = _trace_image(
            flux_orig.shape,
            spline_models,
            spline_scales,
            region_map,
            alpha_orig,
            slope_limit=slope_limit,
            pad=oversample_kwargs["pad"],
        )
        if rotate:
            trace_used = np.rot90(trace_used, k=-1)
            full_trace = np.rot90(full_trace, k=-1)
        model.trace_model = trace_used
        if return_intermediate_models:
            new_models = _intermediate_models(model, [full_trace, trace_used, None, None])
            return model, *new_models
        else:
            return model

    # Oversampled array size
    os_shape = (int(np.ceil(ysize * oversample_factor)), xsize)
    x_os = np.full(os_shape, np.nan)
    y_os = np.full(os_shape, np.nan)

    # Pre-compute coordinates for the new data size
    log.info("Computing oversampled coordinates")
    basex, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))
    newy, oldy = _reindex(0, ysize - 1, scale=oversample_factor)
    y_os[:, :] = oldy[:, None]
    x_os[:, :] = basex[oldy.astype(int), :]
    if mode == "NRS_IFU":
        alpha_os, wave_os = _get_oversampled_coords_nrs_ifu(wcs, x_os, y_os)
    else:
        # Because MIRI was rotated the indexing in the non-rotated frame,
        # the input coordinates need to be adjusted slightly
        det2ab = model.meta.wcs.get_transform("detector", "alpha_beta")
        alpha_os, _, wave_os = det2ab(ysize - y_os - 1, x_os)

    log.info("Oversampling the flux array from the fit trace model")
    flux_os, trace_used, full_trace, linear, residual = oversample_flux(
        flux_orig,
        alpha_orig,
        region_map,
        spline_models,
        spline_scales,
        oversample_factor,
        alpha_os,
        slope_limit=slope_limit,
        psf_optimal=psf_optimal,
        require_ngood=fit_kwargs["require_ngood"],
        **oversample_kwargs,
    )

    log.info("Oversampling error and DQ arrays")
    error_extensions = ["err", "var_rnoise", "var_poisson", "var_flat"]
    errors = {}
    for extname in error_extensions:
        if model.hasattr(extname):
            errors[extname] = getattr(model, extname)
            if rotate:
                errors[extname] = np.rot90(errors[extname])
    dq = model.dq
    if rotate:
        dq = np.rot90(dq)

    # Nearest pixel interpolation for the dq and regions array
    closest_pix = (np.round(y_os).astype(int), np.round(x_os).astype(int))
    dq_os = dq[*closest_pix]
    regions_os = region_map[*closest_pix]

    # Update the DQ image for pixels that used to be NaN, now replaced by spline interpolation.
    # Remove the DO_NOT_USE flag, add FLUX_ESTIMATED
    is_estimated = ~np.isnan(flux_os) & ((dq_os & dqflags.pixel["DO_NOT_USE"]) > 0)
    dq_os[is_estimated] ^= dqflags.pixel["DO_NOT_USE"]
    dq_os[is_estimated] |= dqflags.pixel["FLUX_ESTIMATED"]

    # Simple linear oversample for the error arrays
    errors_os = {}
    for extname, error_array in errors.items():
        error_os = linear_oversample(
            error_array,
            region_map,
            oversample_factor,
            fit_kwargs["require_ngood"],
            edge_limit=0,
            preserve_nan=False,
        )

        # Restore NaNs from the input, except at the estimated locations
        is_nan = ~np.isfinite(error_array[closest_pix])
        error_os[is_nan & ~is_estimated] = np.nan

        # Inflate the errors to account for oversampling covariance
        _inflate_error(error_os, extname, oversample_factor)

        errors_os[extname] = error_os

    # Update the wcs for new pixel scale
    scale_and_shift = Scale(1 / oversample_factor) | Shift(
        -(oversample_factor - 1) / (oversample_factor * 2)
    )
    if mode == "NRS_IFU":
        map_pixels = Identity(1) & scale_and_shift
        model.meta.wcs = _update_wcs_nrs_ifu(model.meta.wcs, map_pixels)
    else:
        # MIRI
        map_pixels = scale_and_shift & Identity(1)
        model.meta.wcs = _update_wcs(model.meta.wcs, map_pixels)

    # If needed, undo all of our rotations before passing back the arrays
    if rotate:
        flux_os = np.rot90(flux_os, k=-1)
        dq_os = np.rot90(dq_os, k=-1)
        wave_os = np.rot90(wave_os, k=-1)
        regions_os = np.rot90(regions_os, k=-1)
        trace_used = np.rot90(trace_used, k=-1)
        full_trace = np.rot90(full_trace, k=-1)
        linear = np.rot90(linear, k=-1)
        residual = np.rot90(residual, k=-1)
        for extname, error_array in errors_os.items():
            errors_os[extname] = np.rot90(error_array, k=-1)

    # Update the model with the oversampled arrays
    model.data = flux_os
    model.dq = dq_os
    model.wavelength = wave_os
    model.trace_model = trace_used
    for extname, error_array in errors_os.items():
        setattr(model, extname, error_array)
    if isinstance(model, datamodels.IFUImageModel):
        model.regions = regions_os

    # Remove some extra arrays if present: no longer needed
    extras = [
        "area",
        "pathloss_point",
        "pathloss_uniform",
        "zeroframe",
    ]
    for name in extras:
        if model.hasattr(name):
            delattr(model, name)

    # Make sure NaNs and DO_NOT_USE flags match in all extensions
    match_nans_and_flags(model)

    # Return intermediate models if needed
    if return_intermediate_models:
        new_models = _intermediate_models(model, [full_trace, trace_used, linear, residual])
        return model, *new_models
    else:
        return model
