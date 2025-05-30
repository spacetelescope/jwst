"""
Utilities for the ATOCA (Darveau-Bernier 2021, in prep).

ATOCA: Algorithm to Treat Order ContAmination (English)
       Algorithme de Traitement d’Ordres ContAmines (French)

@authors: Antoine Darveau-Bernier, Geert Jan Talens
"""

import numpy as np
from numpy.polynomial import Polynomial
import warnings
from scipy.sparse import diags, csr_matrix
from scipy.sparse.linalg import spsolve, lsqr, MatrixRankWarning
from scipy.interpolate import interp1d, RectBivariateSpline, Akima1DInterpolator
from scipy.optimize import minimize_scalar, brentq
from scipy.interpolate import make_interp_spline
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def arange_2d(starts, stops):
    """
    Generate indices on the oversampled wavelength grid.

    Creates a 2D array containing a series of ranges.
    The ranges do not have to be of equal length.

    Parameters
    ----------
    starts : array[int]
        Start values for each range.
    stops : array[int]
        End values for each range.

    Returns
    -------
    out : array[uint16]
        2D array of ranges with invalid values set to -1
    """
    if starts.shape != stops.shape:
        msg = (
            "Shapes of starts and stops are not compatible, "
            "they must either have the same shape or starts must be scalar."
        )
        log.critical(msg)
        raise ValueError(msg)

    if np.any(stops < starts):
        msg = "stops must be everywhere greater or equal to starts."
        log.critical(msg)
        raise ValueError(msg)

    # Compute the length of each range.
    lengths = (stops - starts).astype(int)

    # Initialize the output arrays with invalid value
    nrows = len(stops)
    ncols = np.amax(lengths)
    out = np.ones((nrows, ncols), dtype=np.int16) * -1

    # Compute the indices.
    for irow in range(nrows):
        out[irow, : lengths[irow]] = np.arange(starts[irow], stops[irow])
    return out


def sparse_k(val, k, n_k):
    """
    Transform a 2D array `val` to a sparse matrix.

    Parameters
    ----------
    val : array
        2D array to be transformed
    k : array
        2D array to set column position of values in sparse matrix.
        Negative values used for undefined positions in val.
    n_k : int
        Number of columns in output sparse matrix.

    Returns
    -------
    mat : array
        Sparse matrix to be returned
    """
    # Length of axis 0
    n_i = len(k)

    # Get row index
    i_k = np.indices(k.shape)[0]

    # Take only well defined coefficients
    row = i_k[k >= 0]
    col = k[k >= 0]
    data = val[k >= 0]

    return csr_matrix((data, (row, col)), shape=(n_i, n_k))


def get_wave_p_or_m(wave_map, dispersion_axis=1):
    """
    Compute upper and lower boundaries of a pixel map, given the pixel central value.

    Parameters
    ----------
    wave_map : array[float]
        2d-map of the pixel central wavelength
    dispersion_axis : int, optional
        Which axis is the dispersion axis (0 or 1)

    Returns
    -------
    wave_upper, wave_lower : array[float]
        The wavelength upper and lower boundaries of each pixel, given the central value.
    """
    # Get wavelength boundaries of each pixels
    wave_left, wave_right = _get_wv_map_bounds(wave_map, dispersion_axis=dispersion_axis)

    # The outputs depend on the direction of the spectral axis.
    invalid = wave_map == 0
    if ((wave_right >= wave_left) | invalid).all():
        wave_plus, wave_minus = wave_right, wave_left
    elif ((wave_right <= wave_left) | invalid).all():
        wave_plus, wave_minus = wave_left, wave_right
    else:
        msg = "Some pixels do not follow the expected dispersion axis!"
        log.critical(msg)
        raise ValueError(msg)

    return wave_plus, wave_minus


def _get_wv_map_bounds(wave_map, dispersion_axis=1):
    """
    Compute boundaries of a pixel map, given the pixel central value.

    Parameters
    ----------
    wave_map : array[float]
        2d-map of the pixel central wavelength
    dispersion_axis : int, optional
        Which axis is the dispersion axis (0 or 1)

    Returns
    -------
    wave_top : array[float]
        Wavelength of top edge for each pixel
    wave_bottom : array[float]
        Wavelength of bottom edge for each pixel

    Notes
    -----
    Handling of invalid pixels may lead to unexpected results as follows:
    Bad pixels are completely ignored when computing pixel-to-pixel differences, so
    wv_map=[2,4,6,NaN,NaN,12,14,16] will give wave_top=[1,3,5,0,0,9,13,15]
    because the difference at index 5 was calculated as 12-(12-6)/2=9,
    i.e., as though index 2 and 5 were next to each other.
    A human (or a smarter linear interpolation) would figure out the slope is 2 and
    determine the value of wave_top[5] should most likely be 11.
    This is found not to matter in practice for the current use cases.
    """
    if dispersion_axis == 1:
        # Simpler to use transpose
        wave_map = wave_map.T
    elif dispersion_axis != 0:
        msg = "Dispersion axis must be 0 or 1!"
        log.critical(msg)
        raise ValueError(msg)

    # Initialize arrays.
    wave_top = np.zeros_like(wave_map)
    wave_bottom = np.zeros_like(wave_map)

    # for loop is needed to compute diff in just one spatial direction
    # while skipping invalid values- not trivial to do with array comprehension even
    # using masked arrays
    n_col = wave_map.shape[1]
    for idx in range(n_col):
        wave_col = wave_map[:, idx]

        # Compute the change in wavelength for valid cols
        idx_valid = np.isfinite(wave_col) & (wave_col >= 0)
        wv_col_valid = wave_col[idx_valid]
        delta_wave = np.diff(wv_col_valid) / 2

        # handle edge effects using a constant-difference rule
        delta_wave_top = np.insert(delta_wave, 0, delta_wave[0])
        delta_wave_bottom = np.append(delta_wave, delta_wave[-1])

        # Compute the wavelength values on the top and bottom edges of each pixel.
        wv_col_top = wv_col_valid - delta_wave_top
        wv_col_bottom = wv_col_valid + delta_wave_bottom

        wave_top[idx_valid, idx] = wv_col_top
        wave_bottom[idx_valid, idx] = wv_col_bottom

    # De-transpose if it was transposed for computation
    if dispersion_axis == 1:
        wave_top, wave_bottom = wave_top.T, wave_bottom.T

    return wave_top, wave_bottom


def oversample_grid(wave_grid, n_os):
    """
    Create an oversampled version of the input 1D wavelength grid.

    Parameters
    ----------
    wave_grid : array[float]
        Wavelength grid to be oversampled.
    n_os : int or array[int]
        Oversampling factor. If it is a scalar, take the same value for each
        interval of the grid. If it is an array, n_os specifies the oversampling
        at each interval of the grid, so len(n_os) = len(wave_grid) - 1.

    Returns
    -------
    wave_grid_os : array[float]
        The oversampled wavelength grid.
    """
    # Convert n_os to an array of size len(wave_grid) - 1.
    n_os = np.asarray(n_os)
    if n_os.ndim == 0:
        n_os = np.repeat(n_os, len(wave_grid) - 1)
    elif len(n_os) != (len(wave_grid) - 1):
        msg = "n_os must be a scalar or an array of size len(wave_grid) - 1."
        log.critical(msg)
        raise ValueError(msg)

    # Compute the oversampled grid.
    intervals = 1 / n_os
    intervals = np.insert(np.repeat(intervals, n_os), 0, 0)
    grid = np.cumsum(intervals)
    wave_grid_os = np.interp(grid, np.arange(wave_grid.size), wave_grid)

    # Take only unique values and sort them.
    return np.unique(wave_grid_os)


def _extrapolate_grid(wave_grid, wave_range, poly_ord=1):
    """
    Extrapolate the 1D wavelength grid to cover a given range of values.

    This is done by fitting the derivative with a polynomial of a given order and using it to
    compute subsequent values at both ends of the grid.

    Parameters
    ----------
    wave_grid : array[float]
        Wavelength grid to be extrapolated.
    wave_range : list[float]
        Wavelength range the new grid should cover.
    poly_ord : int
        Order of the polynomial used to fit the derivative of wave_grid.

    Returns
    -------
    wave_grid_ext : array[float]
        The extrapolated 1D wavelength grid.
    """
    if wave_range[0] >= wave_range[-1]:
        msg = "wave_range must be in order [short, long]."
        log.critical(msg)
        raise ValueError(msg)
    if wave_range[0] > wave_grid.max() or wave_range[-1] < wave_grid.min():
        msg = "wave_range must overlap with wave_grid."
        log.critical(msg)
        raise ValueError(msg)
    if wave_range[0] > wave_grid.min() and wave_range[-1] < wave_grid.max():
        return wave_grid

    # Define delta_wave as a function of wavelength by fitting a polynomial.
    delta_wave = np.diff(wave_grid)
    f_delta = Polynomial.fit(wave_grid[:-1], delta_wave, poly_ord).convert()

    # Set a minimum delta value to avoid running forever
    min_delta = delta_wave.min() / 10

    # Extrapolate out-of-bound values on the left-side of the grid.
    grid_left = []
    if wave_range[0] < wave_grid.min():
        # Initialize extrapolated grid with the first value of input grid.
        # This point gets double-counted in the final grid, but then unique is called.
        grid_left = [
            wave_grid.min(),
        ]

        # Iterate until the end of wave_range is reached.
        while True:
            next_delta = f_delta(grid_left[-1])
            next_val = grid_left[-1] - next_delta

            grid_left.append(next_val)
            if next_val < wave_range[0]:
                break
            if next_delta < min_delta:
                raise RuntimeError("Extrapolation failed to converge.")

        # Sort extrapolated vales (and keep only unique).
        grid_left = np.unique(grid_left)

    # Extrapolate out-of-bound values on the right-side of the grid.
    grid_right = []
    if wave_range[-1] > wave_grid.max():
        # Initialize extrapolated grid with the last value of input grid.
        # This point gets double-counted in the final grid, but then unique is called.
        grid_right = [
            wave_grid.max(),
        ]

        # Iterate until the end of wave_range is reached.
        while True:
            next_delta = f_delta(grid_right[-1])
            next_val = grid_right[-1] + next_delta

            grid_right.append(next_val)
            if next_val > wave_range[-1]:
                break
            if next_delta < min_delta:
                raise RuntimeError("Extrapolation failed to converge.")

        # Sort extrapolated values (and keep only unique)
        grid_right = np.unique(grid_right)

    # Combine the extrapolated sections with the original grid.
    return np.concatenate([grid_left, wave_grid, grid_right])


def grid_from_map(wave_map, trace_profile):
    """
    Define a wavelength grid based on the wave_map and trace_profile.

    Takes the wavelength of each column at the center of mass of the spatial profile.

    Parameters
    ----------
    wave_map : array[float]
        Array of the pixel wavelengths for a given order.
    trace_profile : array[float]
        Array of the spatial profile for a given order.

    Returns
    -------
    grid : array[float]
        Output wavelength grid
    cols : array[int]
        Column indices used.
    """
    # Use only valid values by setting weights to zero
    trace_profile[trace_profile < 0] = 0
    trace_profile[wave_map <= 0] = 0

    # handle case where all values are invalid for a given wavelength
    # np.average cannot process sum(weights) = 0, so set them to unity then set NaN afterward
    bad_wls = np.sum(trace_profile, axis=0) == 0
    trace_profile[:, bad_wls] = 1
    center_wv = np.average(wave_map, weights=trace_profile, axis=0)
    center_wv[bad_wls] = np.nan
    center_wv = center_wv[~np.isnan(center_wv)]

    # Make sure the wavelength values are in ascending order.
    return np.sort(center_wv)


def grid_from_map_with_extrapolation(wave_map, trace_profile, wave_range=None, n_os=1):
    """
    Define a wavelength grid based on wave_map, trace_profile, and optional user wave_range.

    Takes the central wavelength at each column given by the center of mass of the spatial profile
    (so one wavelength per column).
    If wave_range is outside of the wave_map, extrapolate.

    Parameters
    ----------
    wave_map : array[float]
        Array of the pixel wavelengths for a given order.
    trace_profile : array[float]
        Array of the spatial profile for a given order.
    wave_range : list[float]
        Minimum and maximum boundary of the grid to generate, in microns.
        Wave_range must include some wavelengths of wave_map.
        Note wave_range is exclusive, in the sense that wave_range[0] and wave_range[1]
        will not be between min(output) and max(output). Instead, min(output) will be
        the smallest value in the extrapolated grid that is greater than wave_range[0]
        and max(output) will be the largest value that is less than wave_range[1].
    n_os : int
        Oversampling of the grid compared to the pixel sampling.

    Returns
    -------
    grid_os : array[float]
        Wavelength grid with oversampling applied
    """
    if wave_map.shape != trace_profile.shape:
        msg = "wave_map and trace_profile must have the same shape."
        log.critical(msg)
        raise ValueError(msg)

    # Different treatment if wave_range is given.
    if wave_range is None:
        grid = grid_from_map(wave_map, trace_profile)
    else:
        # Get an initial estimate of the grid.
        grid = grid_from_map(wave_map, trace_profile)

        # Extrapolate values out of the wv_map if needed
        grid = _extrapolate_grid(grid, wave_range, poly_ord=1)

        # Constrain grid to be within wave_range
        grid = grid[grid >= wave_range[0]]
        grid = grid[grid <= wave_range[-1]]

        # Check if grid and wv_range are compatible
        if len(grid) == 0:
            msg = "Invalid wave_map or wv_range."
            log.critical(msg)
            raise ValueError(msg)

    # Apply oversampling
    return oversample_grid(grid, n_os=n_os)


def _trim_grids(all_grids, grid_range):
    """
    Trim the grids to the wavelength range and remove overlapping parts.

    Remove all parts of the grids that are not in range
    or that are already covered by grids with higher priority,
    i.e., preceding in the list.

    Parameters
    ----------
    all_grids : list[array]
        List of grid (arrays) to trim, in order of importance.
    grid_range : list[float]
        Wavelength range the new grid should cover.

    Returns
    -------
    grids_trimmed : list[array]
        List of trimmed grids.
    """
    grids_trimmed = []
    for grid in all_grids:
        # Remove parts of the grid that are not in the wavelength range
        i_min = np.searchsorted(grid, grid_range[0], side="right")
        i_max = np.searchsorted(grid, grid_range[1], side="left")
        # Make sure it is a valid value and take one grid point past the limit
        # since the oversampling could squeeze some nodes near the limits
        i_min = np.max([i_min - 1, 0])
        i_max = np.min([i_max, len(grid) - 1])
        # Trim the grid
        grid = grid[i_min : i_max + 1]

        # Remove parts of the grid that are already covered
        if len(grids_trimmed) > 0:
            # Use all grids already trimmed (so higher in priority)
            conca_grid = np.concatenate(grids_trimmed)
            # Find values below or above
            is_below = grid < np.min(conca_grid)
            is_above = grid > np.max(conca_grid)

            # Remove values already covered, but keep one
            # index past the limit
            if is_below.any():
                idx = np.max(np.nonzero(is_below))
                idx = np.min([idx + 1, len(grid) - 1])
                grid = grid[: idx + 1]
            if is_above.any():
                idx = np.min(np.nonzero(is_above))
                idx = np.max([idx - 1, 0])
                grid = grid[idx:]

            # If all is covered, no need to do it again, so empty grid.
            if not is_below.any() and not is_above.any():
                grid = np.array([])

        # Save trimmed grid
        grids_trimmed.append(grid)

    return grids_trimmed


def make_combined_adaptive_grid(
    all_grids, all_estimates, grid_range, max_iter=10, rtol=10e-6, max_total_size=1000000
):
    """
    Build an irregular oversampled grid needed to reach a given precision when integrating.

    The grid is built by subdividing iteratively each intervals that
    did not reach the required precision.
    The precision is computed based on the estimate of the integrals
    using a first order Romberg integration.

    Parameters
    ----------
    all_grids : list[array]
        List of grid (arrays) to pass to _adapt_grid, in order of importance.
    all_estimates : list[callable]
        List of function (callable) to estimate the precision needed to oversample the grid.
        Must match the corresponding `grid` in `all_grid`.
    grid_range : list[float]
        Wavelength range the new grid should cover.
    max_iter : int, optional
        Number of times the intervals can be subdivided. The smallest
        subdivison of the grid if max_iter is reached will then be given
        by delta_grid / 2^max_iter. Needs to be greater than zero.
        Default is 10.
    rtol : float, optional
        The desired relative tolerance. Default is 10e-6, so 10 ppm.
    max_total_size : int, optional
        Maximum size of the output grid. Default is 1 000 000.

    Returns
    -------
    os_grid : 1D array
        Oversampled combined grid which minimizes the integration error based on
        Romberg's method
    """
    # Remove unneeded parts of the grids
    all_grids = _trim_grids(all_grids, grid_range)

    # Save native size of each grids (use later to adjust max_grid_size)
    all_sizes = [len(grid) for grid in all_grids]

    # Iterate over grids to build the combined grid
    combined_grid = np.array([])  # Init with empty array
    for i_grid, grid in enumerate(all_grids):
        # Get the max_grid_size, considering the other grids
        # First, remove length already used
        max_grid_size = max_total_size - combined_grid.size
        # Save some space for next grids (at least the native grid size)
        for i_size, size in enumerate(all_sizes):
            if i_size > i_grid:
                max_grid_size = max_grid_size - size
        # Make sure it is at least the size of the native grid.
        max_grid_size = np.max([max_grid_size, all_sizes[i_grid]])

        # Oversample the grid based on tolerance required
        grid, is_converged = _adapt_grid(
            grid, all_estimates[i_grid], max_grid_size=max_grid_size, max_iter=max_iter, rtol=rtol
        )

        # Update grid sizes
        all_sizes[i_grid] = grid.size

        # Check convergence
        if not is_converged:
            msg = "Precision cannot be guaranteed:"
            if grid.size < max_grid_size:
                msg += (
                    f" smallest subdivision 1/{2**max_iter:2.1e}"
                    f" was reached for grid index = {i_grid}"
                )
            else:
                total_size = np.sum(all_sizes)
                msg += " max grid size of "
                msg += " + ".join([f"{size}" for size in all_sizes])
                msg += f" = {total_size} was reached for grid index = {i_grid}."
            log.warning(msg)

        # Combine grids
        combined_grid = np.concatenate([combined_grid, grid])

    # Sort values (and keep only unique).
    # This is necessary because trim_grids allows lowest index of one grid to
    # equal highest index of another grid.
    return np.unique(combined_grid)


def _romberg_diff(b, c, k):
    """
    Compute the differences for the Romberg quadrature corrections.

    See Forman Acton's "Real Computing Made Real," p 143.

    Parameters
    ----------
    b : float or array[float]
        R(n-1, m-1) of Rombergs method.
    c : float or array[float]
        R(n, m-1) of Rombergs method.
    k : int
        The parameter m of Rombergs method.

    Returns
    -------
    R(n, m) : float or array[float]
        Difference between integral estimates of Rombergs method.
    """
    return (4.0**k * c - b) / (4.0**k - 1.0)


def _difftrap(fct, intervals, numtraps):
    """
    Perform part of the trapezoidal rule to integrate a function.

    Assume that we had called difftrap with all lower powers-of-2 starting with 1.
    Calling difftrap only returns the summation of the new ordinates. It does not
    multiply by the width of the trapezoids. This must be performed by the
    caller.

    Note: This function is based on scipy.integrate.quadrature. Adapted to work
    with multiple intervals.

    Parameters
    ----------
    fct : callable
        Function to be integrated.
    intervals : array[float]
        A 2D array of integration intervals of shape (Nx2) or a
        single interval of shape (2,).
    numtraps : int
        The number of trapezoids used to integrate the interval.
        numtraps must be a power of 2.

    Returns
    -------
    s : float
        The sum of function values at the new trapezoid boundaries
        compared to numtraps = numtraps/2. When numtraps = 1 they
        are divided by two.
    """
    # Convert input intervals to numpy array
    intervals = np.asarray(intervals)

    # If intervals is 1D assume it's a single interval.
    if intervals.ndim == 1:
        intervals = intervals[:, np.newaxis]

    # Check the value of numtraps.
    if numtraps <= 0:
        err_msg = "numtraps must be > 0 in difftrap()."
        log.critical(err_msg)
        raise ValueError(err_msg)

    if numtraps == 1:
        # Return the function evaluations for a single trapezoid.
        # Only points at the edge of the interval need to be halved.
        ordsum = 0.5 * (fct(intervals[0]) + fct(intervals[1]))

    elif numtraps % 2:
        err_msg = "numtraps must be a power of 2 in difftrap()."
        log.critical(err_msg)
        raise ValueError(err_msg)
    else:
        # Number of new points compared to lower 2**N multiple of trapezoids.
        numtosum = numtraps / 2

        # Find coordinates of new points.
        h = (intervals[1] - intervals[0]) / numtosum
        lox = intervals[0] + (h * 0.5)
        points = lox[np.newaxis, :] + (h * np.arange(numtosum)[:, np.newaxis])

        # Evaluate and sum the new points.
        ordsum = np.sum(fct(points), axis=0)

    return ordsum


def _estim_integration_err(grid, fct):
    """
    Estimate integration error on each interval of the grid using 1st order Romberg integration.

    Parameters
    ----------
    grid : 1d array [float]
        Grid for integration. Each sections of this grid are treated
        as separate integrals. So if grid has length N; N-1 integrals are
        tested.
    fct : callable
        Function to be integrated.

    Returns
    -------
    err : array[float]
        Absolute error of each integration, with length = length(grid) - 1
    rel_err : array[float]
        Relative error of each integration, with length = length(grid) - 1
    """
    # Change the 1D grid into a 2D set of intervals.
    intervals = np.array([grid[:-1], grid[1:]])
    intrange = np.diff(grid)

    # Estimate of trapezoidal integration without subdivision.
    numtraps = 1
    ordsum = _difftrap(fct, intervals, numtraps)
    trpz = intrange * ordsum / numtraps

    # Estimate with intervals subdivided in 2
    numtraps = 2
    ordsum += _difftrap(fct, intervals, numtraps)
    trpz_sub = intrange * ordsum / numtraps

    # Compute better estimate of the integral
    # using Romberg R(1, 0)
    romb = _romberg_diff(trpz, trpz_sub, 1)

    # Compute errors
    err = np.abs(romb - trpz)
    non_zero = romb != 0
    rel_err = np.full_like(err, np.inf)
    rel_err[non_zero] = np.abs(err[non_zero] / romb[non_zero])

    return err, rel_err


def _adapt_grid(grid, fct, max_grid_size, max_iter=10, rtol=10e-6, atol=1e-6):
    """
    Build an irregular oversampled grid needed to reach a given precision when integrating.

    The grid is built by subdividing iteratively each intervals that
    did not reach the required precision.
    The precision is computed based on the estimate of the integrals
    using a first order Romberg integration.

    See also `scipy.integrate.quadrature.romberg`.

    Parameters
    ----------
    grid : array, required
        Grid for integration. Each sections of this grid are treated
        as separate integrals. So if grid has length N; N-1 integrals are
        optimized.
    fct : callable, required
        Function to be integrated. Must be a function of `grid`
    max_grid_size : int, required
        Maximum size of the output grid.
    max_iter : int, optional
        Number of times the intervals can be subdivided. The smallest
        subdivison of the grid if max_iter is reached will then be given
        by delta_grid / 2^max_iter. Needs to be greater then zero.
        Default is 10.
    rtol : float, optional
        The desired relative tolerance. Default is 10e-6, so 10 ppm.
    atol : float, optional
        The desired absolute tolerance. Default is 1e-6.

    Returns
    -------
    os_grid : 1D array
        Oversampled grid which minimizes the integration error based on Romberg's method
    convergence_flag : bool
        Whether the estimated tolerance was reach everywhere or not.

    References
    ----------
    [1] 'Romberg's method' https://en.wikipedia.org/wiki/Romberg%27s_method
    """
    # Init some flags
    max_size_reached = grid.size > max_grid_size
    if max_size_reached:
        raise ValueError("max_grid_size is too small for the input grid.")

    # Iterate until precision is reached or max_iter
    for _ in range(max_iter):
        # Estimate error using Romberg integration
        abs_err, rel_err = _estim_integration_err(grid, fct)

        # Check where precision is reached
        converged = (rel_err < rtol) | (abs_err < atol)
        is_converged = converged.all()

        # Stop iterating if max grid size was reached
        if max_size_reached or is_converged:
            break

        # Intervals that didn't reach the precision will be subdivided
        n_oversample = np.full(rel_err.shape, 2, dtype=int)
        # No subdivision for the converged ones
        n_oversample[converged] = 1

        # Check if the maximum size will be reached.
        # If so, prioritize the intervals with the largest estimated errors
        # to reach the maximum size
        os_grid_size = n_oversample.sum()
        if os_grid_size > max_grid_size:
            max_size_reached = True

            # How many nodes can be added to reach max?
            n_nodes_remaining = max_grid_size - grid.size

            # Find the position of the nodes with the largest error
            idx_largest_err = np.argsort(rel_err)[-n_nodes_remaining:]

            # Build new oversample array and assign only largest errors
            n_oversample = np.ones(rel_err.shape, dtype=int)
            n_oversample[idx_largest_err] = 2

        # Generate oversampled grid (subdivide). Returns sorted and unique grid.
        grid = oversample_grid(grid, n_os=n_oversample)

    return grid, is_converged


def throughput_soss(wavelength, throughput):
    """
    Create an interpolator for the throughput values.

    Parameters
    ----------
    wavelength : array[float]
        A wavelength array.
    throughput : array[float]
        The throughput values corresponding to the wavelengths.

    Returns
    -------
    interpolator : callable
        A function that interpolates the throughput values. Accepts an array
        of wavelengths and returns the interpolated throughput values.

    Notes
    -----
    Throughput is always zero at min, max of wavelength.
    """
    wavelength = np.sort(wavelength)
    wl_min, wl_max = np.min(wavelength), np.max(wavelength)
    throughput[0] = 0.0
    throughput[-1] = 0.0
    interp = make_interp_spline(wavelength, throughput, k=3, bc_type=("clamped", "clamped"))

    def interpolator(wv):
        wv = np.clip(wv, wl_min, wl_max)
        return interp(wv)

    return interpolator


class WebbKernel:
    """The JWST kernel."""

    def __init__(self, wave_kernels, kernels, wave_trace, n_pix):
        """
        Initialize the kernel object.

        Parameters
        ----------
        wave_kernels : array[float]
            Wavelength array for the kernel. Must have same shape as kernels.
        kernels : array[float]
            Kernel for throughput array.
            Dimensions are (wavelength, oversampled pixels).
            Center (~max throughput) of the kernel is at the center of the 2nd axis.
        wave_trace : array[float]
            1-D trace of the detector central wavelengths for the given order.
            Since kernels are originally defined in the pixel space, this is used to
            convert to wavelength space.
        n_pix : int
            Number of detector pixels spanned by the kernel. Second axis of kernels
            has shape (n_os * n_pix) - (n_os - 1), where n_os is the
            spectral oversampling factor.
        """
        self.n_pix = n_pix

        # Mask where trace is equal to 0
        wave_trace = np.ma.array(wave_trace, mask=(wave_trace == 0))

        # Force trace to have the red wavelengths at the end of the detector
        if np.diff(wave_trace).mean() < 0:
            wave_trace = np.flip(wave_trace)

        # Create oversampled pixel position array. Center index should have value 0.
        self.pixels = np.linspace(-(n_pix // 2), n_pix // 2, wave_kernels.shape[0])

        # `wave_kernel` has only the value of the central wavelength
        # of the kernel at each points because it's a function
        # of the pixels (so depends on wv solution).
        wave_center = wave_kernels[0, :]

        # Use the wavelength solution to create a mapping between pixels and wavelengths
        wave_min = np.amin(wave_trace[wave_trace > 0])
        wave_max = np.amax(wave_trace[wave_trace > 0])
        i_min = np.searchsorted(wave_center, wave_min)
        i_max = np.searchsorted(wave_center, wave_max) - 1

        # i_min, i_max correspond to the min, max indices of the kernel that are represented
        # on the detector. Use those to define the boundaries of the interpolation to use
        # in the RectBivariateSpline interpolation
        bbox = [
            None,
            None,
            wave_center[np.maximum(i_min - 1, 0)],
            wave_center[np.minimum(i_max + 1, len(wave_center) - 1)],
        ]

        # Keep only kernels that fall on the detector.
        self.kernels = kernels[:, i_min : i_max + 1].copy()
        wave_kernels = wave_kernels[:, i_min : i_max + 1].copy()
        wave_center = np.array(wave_kernels[0])

        # Save minimum kernel value (greater than zero)
        self.min_value = np.min(self.kernels[(self.kernels > 0.0)])

        # Then find the pixel closest to each kernel center
        # and use the surrounding pixels (columns)
        # to get the wavelength. At the boundaries,
        # wavelength might not be defined or falls out of
        # the detector, so fit a 1-order polynomial to
        # extrapolate. The polynomial is also used to interpolate
        # for oversampling.
        i_surround = np.arange(-(n_pix // 2), n_pix // 2 + 1)
        poly = []
        for i_cen, wv_c in enumerate(wave_center):
            wv = np.ma.masked_all(i_surround.shape)

            # Closest pixel wv
            i_col = np.argmin(np.abs(wave_trace - wv_c))
            # Update wavelength center value
            # (take the nearest pixel center value)
            wave_center[i_cen] = wave_trace[i_col]

            # Surrounding columns
            index = i_col + i_surround

            # Make sure it's on the detector
            i_good = (index >= 0) & (index < wave_trace.size)

            # Assign wv values
            wv[i_good] = wave_trace[index[i_good]]

            # Fit n=1 polynomial
            f = Polynomial.fit(i_surround[~wv.mask], wv[~wv.mask], 1).convert()
            poly_i = f.coef[::-1]  # Reverse order to match old behavior from legacy np.polyval

            # Project on os pixel grid
            wave_kernels[:, i_cen] = f(self.pixels)

            # Save coeffs
            poly.append(poly_i)

        # Save computed attributes
        self.wave_kernels = wave_kernels
        self.wave_center = wave_center
        self.poly = np.array(poly)

        # 2D Interpolate
        self.f_ker = RectBivariateSpline(self.pixels, self.wave_center, self.kernels, bbox=bbox)

    def __call__(self, wave, wave_c):
        """
        Return the kernel value, given the wavelength and the kernel central wavelength.

        Wavelengths that are out of bounds will be extrapolated.

        Parameters
        ----------
        wave : array[float]
            Wavelength where the kernel is projected.
        wave_c : array[float]
            Central wavelength of the kernel.

        Returns
        -------
        out : array[float]
            The kernel value.
        """
        wave_center = self.wave_center
        poly = self.poly
        n_wv_c = len(wave_center)

        # First, convert wavelength value into pixels using self.poly to interpolate

        # Find corresponding interval
        i_wv_c = np.searchsorted(wave_center, wave_c) - 1

        # Extrapolate values out of bounds
        i_wv_c[i_wv_c < 0] = 0
        i_wv_c[i_wv_c >= (n_wv_c - 1)] = n_wv_c - 2

        # Compute coefficients that interpolate along wv_centers
        d_wv_c = wave_center[i_wv_c + 1] - wave_center[i_wv_c]
        a_c = (wave_center[i_wv_c + 1] - wave_c) / d_wv_c
        b_c = (wave_c - wave_center[i_wv_c]) / d_wv_c

        # Compute a_pix and b_pix from the equation:
        # pix = a_pix * lambda + b_pix
        a_pix = 1 / (a_c * poly[i_wv_c, 0] + b_c * poly[i_wv_c + 1, 0])
        b_pix = -(a_c * poly[i_wv_c, 1] + b_c * poly[i_wv_c + 1, 1])
        b_pix /= a_c * poly[i_wv_c, 0] + b_c * poly[i_wv_c + 1, 0]

        # Compute pixel values
        pix = a_pix * wave + b_pix

        # Second, compute kernel value on the interpolation grid (pixel x wv_center)
        webbker = self.f_ker(pix, wave_c, grid=False)

        # Make sure it's not negative and greater than the min value,
        # set pixels outside range to zero
        webbker = np.clip(webbker, self.min_value, None)
        webbker[pix > self.n_pix // 2] = 0
        webbker[pix < -(self.n_pix // 2)] = 0

        return webbker


def _constant_kernel_to_2d(c, grid_range):
    """
    Build a 2D kernel array with a constant 1D kernel as input.

    Parameters
    ----------
    c : float or size-1 np.ndarray
        Constant value to expand into a 2-D kernel
    grid_range : list[int]
        Indices over which convolution is defined on grid.

    Returns
    -------
    kernel_2d : array[float]
        2D array of input 1D kernel tiled over axis with
        length equal to difference of grid_range values.
    """
    # Assign range where the convolution is defined on the grid
    a, b = grid_range

    # Get length of the convolved axis
    n_k_c = b - a

    # Return a 2D array with this length
    return np.tile(np.atleast_1d(c), (n_k_c, 1)).T


def _get_wings(fct, grid, h_len, i_a, i_b):
    """
    Compute values of the kernel at grid[+-h_len].

    Parameters
    ----------
    fct : callable
        Function that returns the value of the kernel, given
        a grid value and the center of the kernel.
        fct(grid, center) = kernel
        grid and center have the same length.
    grid : array[float]
        Grid where the kernel is projected
    h_len : int
        Half-length where we compute kernel value.
    i_a : int
        Index of grid axis 0 where to apply convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[i_a:i_b].
    i_b : int
        Index of grid axis 1 where to apply convolution.

    Returns
    -------
    left : array[float]
        Kernel values at left wing.
    right : array[float]
        Kernel values at right wing.
    """
    # Save length of the non-convolved grid
    n_k = len(grid)

    # Get length of the convolved axis
    n_k_c = i_b - i_a

    # Init values
    left, right = np.zeros(n_k_c), np.zeros(n_k_c)

    # Add the left value on the grid
    # Possibility that it falls out of the grid;
    # take first value of the grid if so.
    i_grid = np.max([0, i_a - h_len])

    # Save the new grid
    grid_new = grid[i_grid : i_b - h_len]

    # reuse dummy variable `i_grid`
    i_grid = len(grid_new)

    # Compute kernel at the left end.
    # `i_grid` accounts for smaller length.
    ker = fct(grid_new, grid[i_b - i_grid : i_b])
    left[-i_grid:] = ker

    # Add the right value on the grid
    # Possibility that it falls out of the grid;
    # take last value of the grid if so.
    # Same steps as the left end (see above)
    i_grid = np.min([n_k, i_b + h_len])
    grid_new = grid[i_a + h_len : i_grid]
    i_grid = len(grid_new)
    ker = fct(grid_new, grid[i_a : i_a + i_grid])
    right[:i_grid] = ker

    return left, right


def _trpz_weight(grid, length, shape, i_a, i_b):
    """
    Compute weights due to trapezoidal integration.

    Parameters
    ----------
    grid : array[float]
        Grid where the integration is projected
    length : int
        Length of the kernel
    shape : tuple[int]
        Shape of the compact convolution 2d array
    i_a : int
        Index of grid axis 0 where to apply convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[i_a:i_b].
    i_b : int
        Index of grid axis 1 where to apply convolution.

    Returns
    -------
    out : array[float]
        2D array with shape according to input shape
    """
    # Index of each element on the convolution matrix
    # with respect to the non-convolved grid
    # `i_grid` has the shape (N_k_convolved, kernel_length - 1)
    i_grid = np.indices(shape)[0] - (length // 2)
    i_grid = np.arange(i_a, i_b)[None, :] + i_grid[:-1, :]

    # Set values out of grid to -1
    i_bad = (i_grid < 0) | (i_grid >= len(grid) - 1)
    i_grid[i_bad] = -1

    # Delta lambda
    d_grid = np.diff(grid)

    # Compute weights from trapezoidal integration
    weight = 0.5 * d_grid[i_grid]
    weight[i_bad] = 0

    # Fill output
    out = np.zeros(shape)
    out[:-1] += weight
    out[1:] += weight

    return out


def _fct_to_array(fct, grid, grid_range, thresh):
    """
    Build a compact kernel 2d array based on a kernel function and a grid to project the kernel.

    Parameters
    ----------
    fct : callable
        Function that returns the value of the kernel, given
        a grid value and the center of the kernel.
        fct(grid, center) = kernel
        grid and center have the same length.
    grid : array[float]
        Grid where the kernel is projected
    grid_range : list[int] or tuple[int]
        Indices of the grid where to apply the convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[grid_range[0]:grid_range[1]].
    thresh : float, required
        Threshold to define the maximum length of the kernel.
        Truncate when `kernel` < `thresh`.

    Returns
    -------
    kern_array : array[float]
        2D array of kernel projected onto grid.
    """
    # Assign range where the convolution is defined on the grid
    i_a, i_b = grid_range

    # Init 2-D array with first dimension length 1, with the value at kernel's center
    out = fct(grid, grid)[i_a:i_b][np.newaxis, ...]

    # Add wings: Generate a 2D array of the grid iteratively until
    # thresh is reached everywhere.
    length = 1
    h_len = 0  # Half length
    while True:
        h_len += 1

        # Compute next left and right ends of the kernel
        left, right = _get_wings(fct, grid, h_len, i_a, i_b)

        # Check if they are all below threshold.
        if (left < thresh).all() and (right < thresh).all():
            break  # Stop iteration
        else:
            # Update kernel length
            length += 2

            # Set value to zero if smaller than threshold
            left[left < thresh] = 0.0
            right[right < thresh] = 0.0

            # add new values to output
            out = np.vstack([left, out, right])

    # Weights due to integration (from the convolution)
    weights = _trpz_weight(grid, length, out.shape, i_a, i_b)

    return out * weights


def _sparse_c(ker, n_k, i_zero):
    """
    Convert a convolution kernel in compact form (N_ker, N_k_c) to sparse form (N_k_c, N_k).

    N_k_c represents the length of the convolved grid, N_k the length of the original grid.

    Parameters
    ----------
    ker : array[float]
        Convolution kernel with shape (N_kernel, N_kc)
    n_k : int
        Length of the original grid
    i_zero : int
        Position of the first element of the convolved grid
        in the original grid.

    Returns
    -------
    matrix : array[float]
        Sparse form of the input convolution kernel
    """
    # Assign kernel length and convolved axis length
    n_ker, n_k_c = ker.shape

    # Algorithm works for odd kernel grid
    if n_ker % 2 != 1:
        err_msg = "Length of the convolution kernel given to _sparse_c should be odd."
        log.critical(err_msg)
        raise ValueError(err_msg)

    # Assign half-length
    h_len = (n_ker - 1) // 2

    # Define each diagonal of the sparse convolution matrix
    diag_val, offset = [], []
    for i_ker, i_k_c in enumerate(range(-h_len, h_len + 1)):
        i_k = i_zero + i_k_c

        if i_k < 0:
            diag_val.append(ker[i_ker, -i_k:])
        else:
            diag_val.append(ker[i_ker, :])

        offset.append(i_k)

    # Build convolution matrix
    return diags(diag_val, offset, shape=(n_k_c, n_k), format="csr")


def get_c_matrix(kernel, grid, i_bounds=None, thresh=1e-5):
    """
    Return a convolution matrix.

    Returns a sparse matrix (N_k_convolved, N_k).
    N_k is the length of the grid on which the convolution
    will be applied, N_k_convolved is the length of the
    grid after convolution and N_ker is the maximum length of
    the kernel.
    The convolution can be applied on an array f | f = fct(grid)
    by a simple matrix multiplication:
    f_convolved = c_matrix.dot(f)

    Parameters
    ----------
    kernel : ndarray (2D) or callable
        Convolution kernel. Can be already 2D (N_ker, N_k_convolved),
        giving the kernel for each items of the convolved grid.
        Can be a callable
        with the form f(x, x0) where x0 is the position of the center of
        the kernel. Must return a 1D array (len(x)), so a kernel value
        for each pairs of (x, x0).
    grid : 1D np.array
        The grid on which the convolution will be applied.
        For example, if C is the convolution matrix,
        f_convolved = C.f(grid)
    i_bounds : 2-elements object, optional, default None
        The bounds of the grid on which the convolution is defined.
        For example, if bounds = (a,b),
        then grid_convolved = grid[a <= grid <= b].
        It dictates also the dimension of f_convolved.
        If None, the convolution is defined on the whole grid.
    thresh : float, optional
        Only used when `kernel` is callable to define the maximum
        length of the kernel. Truncate when `kernel` < `thresh`

    Returns
    -------
    c_matrix : array[float]
        Convolution matrix in sparse form (N_k_convolved, N_k).
    """
    # Define range where the convolution is defined on the grid.
    if i_bounds is None:
        a, b = 0, len(grid)

    else:
        # Make sure it is absolute index, not relative
        # So no negative index.
        if i_bounds[1] < 0:
            i_bounds[1] = len(grid) + i_bounds[1]

        a, b = i_bounds

    # Generate a 2D kernel of shape (N_kernel x N_kc)
    if callable(kernel):
        kernel = _fct_to_array(kernel, grid, [a, b], thresh)

    elif kernel.size == 1:
        kernel = _constant_kernel_to_2d(kernel, [a, b])

    elif kernel.ndim != 2:
        msg = "Input kernel to get_c_matrix must be callable or2-dimensional array."
        log.critical(msg)
        raise ValueError(msg)

    # Normalize
    kernel = kernel / np.nansum(kernel, axis=0)

    # Convert to a sparse matrix.
    return _sparse_c(kernel, len(grid), a)


def _finite_diff(x):
    """
    Return the finite difference matrix operator based on x.

    Parameters
    ----------
    x : array[float]
        Input array

    Returns
    -------
    diff_matrix : array[float]
        Sparse matrix. When applied to x `diff_matrix.dot(x)`,
        the result is the same as np.diff(x)
    """
    n_x = len(x)
    diff_matrix = diags([-1.0], shape=(n_x - 1, n_x))
    diff_matrix += diags([1.0], 1, shape=(n_x - 1, n_x))
    return diff_matrix


def finite_first_d(grid):
    """
    Return the first derivative operator based on grid.

    Parameters
    ----------
    grid : array[float]
        Grid where the first derivative will be computed.

    Returns
    -------
    first_d : array[float]
        Operator to compute the first derivative, so that
        f' = first_d.dot(f), where f is a function
        projected on `grid`.
    """
    # Finite difference operator
    d_matrix = _finite_diff(grid)

    # Delta lambda
    d_grid = d_matrix.dot(grid)

    # First derivative operator
    return diags(1.0 / d_grid).dot(d_matrix)


def _curvature_finite(factors, log_reg2, log_chi2):
    """
    Compute the curvature in log space using finite differences.

    Parameters
    ----------
    factors : array[float]
        Regularisation factors (not in log).
    log_reg2 : array[float]
        Norm-2 of the regularisation term (in log10).
    log_chi2 : array[float]
        Norm-2 of the chi2 term (in log10).

    Returns
    -------
    factors : array[float]
        Sorted and cut version of input factors array.
    curvature : array[float]
        Second derivative of the log10 of the regularized chi2
    """
    # Make sure it is sorted according to the factors
    idx = np.argsort(factors)
    factors, log_chi2, log_reg2 = factors[idx], log_chi2[idx], log_reg2[idx]

    # Get first and second derivatives
    chi2_deriv = _get_finite_derivatives(factors, log_chi2)
    reg2_deriv = _get_finite_derivatives(factors, log_reg2)

    # Compute the curvature according to Hansen 2001
    #
    # Numerator of the curvature
    numerator = chi2_deriv[0] * reg2_deriv[1]
    numerator -= reg2_deriv[0] * chi2_deriv[1]
    # Denominator of the curvature
    denom = reg2_deriv[0] ** 2 + chi2_deriv[0] ** 2
    # Combined
    curv = 2 * numerator / np.power(denom, 1.5)

    # Since the curvature is not define at the ends of the array,
    # cut the factors array
    factors = factors[1:-1]

    return factors, curv


def _get_finite_derivatives(x_array, y_array):
    """
    Compute first and second finite derivatives.

    Parameters
    ----------
    x_array : array[float]
        1D array of x values.
    y_array : array[float]
        1D array of y values.

    Returns
    -------
    mean_first_d : array[float]
        Mean of left and right finite derivatives
    second_d : array[float]
        Second finite derivative
    """
    # Compute first finite derivative
    first_d = np.diff(y_array) / np.diff(x_array)
    # Take the mean of the left and right derivative
    mean_first_d = 0.5 * (first_d[1:] + first_d[:-1])

    # Compute second finite derivative
    second_d = 0.5 * np.diff(first_d) / (x_array[2:] - x_array[:-2])

    return mean_first_d, second_d


def _get_interp_idx_array(idx, relative_range, max_length):
    """
    Generate array given the relative range around an index.

    Parameters
    ----------
    idx : int
        Center index value
    relative_range : iterable[int]
        Relative bounds around center value to create new array
    max_length : int
        Upper bound on range of indices to provide

    Returns
    -------
    array[int]
        Output array of indices
    """
    # Convert to absolute index range
    abs_range = [idx + d_idx for d_idx in relative_range]

    # Make sure it's still a valid index
    abs_range[0] = np.max([abs_range[0], 0])
    abs_range[-1] = np.min([abs_range[-1], max_length])

    # Convert to slice
    return np.arange(*abs_range, 1)


def _minimize_on_grid(factors, val_to_minimize, interpolate=True, interp_index=None):
    """
    Find minimum of a grid using akima spline interpolation to get a finer estimate.

    Parameters
    ----------
    factors : array[float]
        1D array of Tikhonov factors for which value array is calculated
    val_to_minimize : array[float]
        1D array of values to be minimized, e.g. chi^2 or curvature.
    interpolate : bool, optional
        If True, use akima spline interpolation to find a finer minimum;
        otherwise, return minimum value in array. Default is true.
    interp_index : iterable[int], optional
        Relative range of grid indices around the minimum value to interpolate
        across. If not specified, defaults to [-2,4].

    Returns
    -------
    min_fac : float
        The factor with minimized error/curvature.
    """
    if interp_index is None:
        interp_index = [-2, 4]

    # Only keep finite values
    idx_finite = np.isfinite(val_to_minimize)
    factors = factors[idx_finite]
    val_to_minimize = val_to_minimize[idx_finite]

    # Get position the minimum
    idx_min = np.argmin(val_to_minimize)

    # If the min is on the one of the boundary, then do not interpolate
    if idx_min == 0 or idx_min == (len(val_to_minimize) - 1):
        interpolate = False

    if interpolate:
        # Interpolate to get a finer estimate
        # Une index only around the best value
        max_length = len(val_to_minimize)
        index = _get_interp_idx_array(idx_min, interp_index, max_length)

        # Akima spline in log space
        x_val, y_val = np.log10(factors[index]), val_to_minimize[index]
        i_sort = np.argsort(x_val)
        x_val, y_val = x_val[i_sort], y_val[i_sort]
        fct = Akima1DInterpolator(x_val, y_val)

        # Find min
        bounds = (x_val.min(), x_val.max())
        opt_args = {"bounds": bounds, "method": "bounded"}
        min_fac = minimize_scalar(fct, **opt_args).x

        # Back to linear scale
        min_fac = 10.0**min_fac

    else:
        # Simply return the min value
        # if no interpolation required
        min_fac = factors[idx_min]

    return min_fac


def _find_intersect(factors, y_val, thresh, interpolate=True, search_range=None):
    """
    Find the root of y_val - thresh (so the intersection between thresh and y_val).

    Parameters
    ----------
    factors : array[float]
        1D array of Tikhonov factors for which value array is calculated
    y_val : array[float]
        1D array of values.
    thresh : float
        Threshold use in 'd_chi2' mode. Find the highest factor where the
        derivative of the chi2 derivative is below thresh.
    interpolate : bool, optional, default True
        If True, use interpolation to find a finer minimum;
        otherwise, return minimum value in array.
    search_range : iterable[int], optional, default [0,3]
        Relative range of grid indices around the value to interpolate.

    Returns
    -------
    float
        Factor corresponding to the best approximation of the intersection point.
    """
    if search_range is None:
        search_range = [0, 3]

    # Only keep finite values
    idx_finite = np.isfinite(y_val)
    factors = factors[idx_finite]
    y_val = y_val[idx_finite]

    # Make sure sorted
    idx_sort = np.argsort(factors)
    factors, y_val = factors[idx_sort], y_val[idx_sort]

    # Check if the threshold is reached
    cond_below = y_val < thresh
    if cond_below.any():
        # Find where the threshold is crossed
        idx_below = np.where(cond_below)[0]
        # Take the largest index (so the highest factor)
        idx_below = np.max(idx_below)
        # If it happens to be the last element of the array...
        if idx_below == (len(factors) - 1):
            # ... no need to interpolate
            interpolate = False
    else:
        # Take the lowest factor value
        idx_below = 0
        # No interpolation needed
        interpolate = False

    if interpolate:
        # Interpolate with log10(factors) to get a finer estimate
        x_val = np.log10(factors)
        d_chi2_spl = interp1d(x_val, y_val - thresh, kind="linear")

        # Use index only around the best value
        max_length = len(y_val)
        index = _get_interp_idx_array(idx_below, search_range, max_length)

        # Find the root
        bracket = (x_val[index[0]], x_val[index[-1]])
        best_val = brentq(d_chi2_spl, *bracket)

        # Back to linear scale
        best_val = 10.0**best_val

    else:
        # Simply return the value
        best_val = factors[idx_below]

    return best_val


def _soft_l1(z):
    return 2 * ((1 + z) ** 0.5 - 1)


def _cauchy(z):
    return np.log(1 + z)


def _linear(z):
    return z


LOSS_FUNCTIONS = {"soft_l1": _soft_l1, "cauchy": _cauchy, "linear": _linear}
DEFAULT_THRESH_DERIVATIVE = {"chi2": 1e-5, "chi2_soft_l1": 1e-4, "chi2_cauchy": 1e-3}


class TikhoTests(dict):
    """Store results of Tikhonov tests for different factors."""

    def __init__(self, test_dict, default_chi2="chi2_cauchy"):
        """
        Merge output of Tikhonov solver with chi2 and curvature.

        Parameters
        ----------
        test_dict : dict
            Dictionary holding arrays for `factors`, `solution`, `error`, `reg`, and `grid`.
        default_chi2 : str
            Type of chi2 loss used by default. Options are chi2, chi2_soft_l1, chi2_cauchy.
        """
        # Define the number of data points
        # (length of the "b" vector in the tikhonov regularisation)
        n_points = len(test_dict["error"][0].squeeze())

        # Save attributes
        self.n_points = n_points
        self.default_chi2 = default_chi2
        self.default_thresh = DEFAULT_THRESH_DERIVATIVE

        # Initialize so it behaves like a dictionary
        super().__init__(test_dict)

        chi2_loss = {"chi2": "linear", "chi2_soft_l1": "soft_l1", "chi2_cauchy": "cauchy"}
        for chi2_type, loss in chi2_loss.items():
            try:
                # Save the chi2
                self[chi2_type]
            except KeyError:
                self[chi2_type] = self._compute_chi2(loss)

    def _compute_chi2(self, loss):
        """
        Calculate the reduced chi squared statistic.

        Parameters
        ----------
        loss : str
            Type of loss function to use. Options are 'linear', 'soft_l1', 'cauchy'.

        Returns
        -------
        float
            Sum of the squared error array divided by the number of data points
        """
        try:
            loss = LOSS_FUNCTIONS[loss]
        except KeyError:
            msg = (
                f"loss={loss} not a valid key. "
                f"Must be one of {[LOSS_FUNCTIONS.keys()]} or callable."
            )
            raise KeyError(msg) from None

        # Compute the reduced chi^2 for all tests
        chi2 = np.nanmean(loss(self["error"] ** 2), axis=-1)
        # Remove residual dimensions
        return chi2.squeeze()

    def _get_chi2_derivative(self):
        """
        Compute derivative of the chi2 with respect to log10(factors).

        Returns
        -------
        factors_leftd : array[float]
            Factors array, shortened to match length of derivative.
        d_chi2 : array[float]
            Derivative of chi squared array with respect to log10(factors)
        """
        key = self.default_chi2

        # Compute finite derivative
        fac_log = np.log10(self["factors"])
        d_chi2 = np.diff(self[key]) / np.diff(fac_log)

        # Update size of factors to fit derivatives
        # Equivalent to derivative on the left side of the nodes
        factors_leftd = self["factors"][1:]

        return factors_leftd, d_chi2

    def _compute_curvature(self):
        """
        Compute the curvature of the l-curve in log-log space.

        Returns
        -------
        factors : array[float]
            Regularisation factors
        curvature : array[float]
            Curvature of the l-curve
        """
        key = self.default_chi2

        # Compute the curvature...
        # Get the norm-2 of the regularisation term
        reg2 = np.nansum(self["reg"] ** 2, axis=-1)

        return _curvature_finite(self["factors"], np.log10(self[key]), np.log10(reg2))

    def best_factor(self, mode="curvature"):
        """
        Compute the best scale factor for Tikhonov regularisation.

        Best factor is determined by taking the factor giving the highest logL on
        the detector or the highest curvature of the l-curve,
        depending on the chosen mode.

        Parameters
        ----------
        mode : str
            How to find the best factor: 'chi2', 'curvature' or 'd_chi2'.

        Returns
        -------
        float
            Best scale factor as determined by the selected algorithm
        """
        key = self.default_chi2
        thresh = self.default_thresh[key]

        # Number of factors
        n_fac = len(self["factors"])

        # Determine the mode (what do we minimize?)
        if mode == "curvature" and n_fac > 2:
            # Compute the curvature
            factors, curv = self._compute_curvature()

            # Find min factor
            best_fac = _minimize_on_grid(factors, curv)

        elif mode == "chi2":
            # Simply take the chi2 and factors
            factors = self["factors"]
            y_val = self[key]

            # Find min factor
            best_fac = _minimize_on_grid(factors, y_val)

        elif mode == "d_chi2" and n_fac > 1:
            # Compute the derivative of the chi2
            factors, y_val = self._get_chi2_derivative()

            # Remove values for the higher factors that
            # are not already below thresh. If not _find_intersect
            # would just return the last value of factors.
            i_last = -1
            while abs(i_last) <= len(y_val):
                # Check derivative
                if y_val[i_last] > thresh:
                    # Save index in slice
                    idx = slice(0, i_last)
                    # break so `else` will be skipped
                    break
                # Update index
                i_last -= 1

            # If all the values were passed without breaking,
            # do not remove any values
            else:
                idx = slice(None)

            # Find intersection with threshold
            best_fac = _find_intersect(factors[idx], y_val[idx], thresh)

        elif mode in ["curvature", "d_chi2", "chi2"]:
            best_fac = np.max(self["factors"])
            msg = (
                f"Could not compute {mode} because number of factors {n_fac} "
                "is too small for that mode."
                f"Setting best factor to max factor: {best_fac:.5e}"
            )
            log.warning(msg)

        else:
            msg = f"`mode`={mode} is not a valid option for TikhoTests.best_factor()."
            log.critical(msg)
            raise ValueError(msg)

        # Return estimated best scale factor
        return best_fac


def try_solve_two_methods(matrix, result):
    """
    Solve sparse matrix equation A.x=b, reverting to least-squared solver when spsolve fails.

    On rare occasions spsolve's approximation of the matrix is not appropriate
    and fails on good input data.

    Parameters
    ----------
    matrix : array-like
        Matrix A in the system to solve A.x = b
    result : array-like
        Vector b in the system to solve A.x = b

    Returns
    -------
    array
        Solution x of the system (1d array)
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(action="error", category=MatrixRankWarning)
        try:
            return spsolve(matrix, result)
        except MatrixRankWarning:
            log.info("ATOCA matrix solve failed with spsolve. Retrying with least-squares.")
            return lsqr(matrix, result)[0]


class Tikhonov:
    """
    Use Tikhonov regularization to solve the ill-posed problem A.x = b.

    The matrix A is accidentally singular or close to singularity. Tikhonov regularization
    adds a regularization term in the equation and aim to minimize the
    equation: ||A.x - b||^2 + ||gamma.x||^2
    where gamma is the Tikhonov regularization matrix.
    """

    def __init__(self, a_mat, b_vec, t_mat):
        """
        Initialize the solver.

        Parameters
        ----------
        a_mat : matrix-like object (2d)
            Matrix A in the system to solve A.x = b
        b_vec : vector-like object (1d)
            Vector b in the system to solve A.x = b
        t_mat : matrix-like object (2d)
            Tikhonov regularisation matrix to be applied on b_vec.
        """
        # Save input matrix
        self.a_mat = a_mat
        self.b_vec = b_vec
        self.t_mat = t_mat

        # Pre-compute some matrix for the linear system to solve
        self.t_mat_2 = (t_mat.T).dot(t_mat)  # squared tikhonov matrix
        self.a_mat_2 = a_mat.T.dot(a_mat)  # squared model matrix
        self.result = (a_mat.T).dot(b_vec.T)
        self.idx_valid = (self.result.toarray() != 0).squeeze()  # valid indices to use

        # Save other attributes
        self.test = None

    def solve(self, factor=1.0):
        """
        Solve the Tikhonov regularization problem.

        Minimize the equation ||A.x - b||^2 + ||gamma.x||^2
        by solving (A_T.A + gamma_T.gamma).x = A_T.b
        gamma is the Tikhonov matrix multiplied by a scale factor

        Parameters
        ----------
        factor : float, optional
            Multiplicative constant of the regularization matrix

        Returns
        -------
        array
            Solution of the system (1d array)
        """
        # Get needed attributes
        a_mat_2 = self.a_mat_2
        result = self.result
        t_mat_2 = self.t_mat_2
        idx = self.idx_valid

        # Matrix gamma squared (with scale factor)
        gamma_2 = factor**2 * t_mat_2

        # Finalize building matrix
        matrix = a_mat_2 + gamma_2

        # Initialize solution
        solution = np.full(matrix.shape[0], np.nan)

        # Solve
        matrix = matrix[idx, :][:, idx]
        result = result[idx]
        solution[idx] = try_solve_two_methods(matrix, result)

        return solution

    def test_factors(self, factors):
        """
        Test multiple candidate Tikhonov factors.

        Parameters
        ----------
        factors : array[float]
            1D array of factors to test

        Returns
        -------
        TikhoTests
            Dictionary of test results
        """
        log.info("Testing factors...")

        # Get relevant attributes
        b_vec = self.b_vec
        a_mat = self.a_mat
        t_mat = self.t_mat

        # Init outputs
        sln, err, reg = [], [], []

        # Test all factors
        for i_fac, factor in enumerate(factors):
            # Save solution
            sln.append(self.solve(factor))

            # Save error A.x - b
            this_err = a_mat.dot(sln[-1]) - b_vec
            # initially this is a np.matrix of shape (1, n_pixels); flatten and make array
            err.append(np.array(this_err).flatten())

            # Save regularization term
            reg_i = t_mat.dot(sln[-1])
            reg.append(reg_i)

            # Print
            message = f"{i_fac}/{len(factors)}"
            log.info(message)

        # Final message output
        message = f"{i_fac + 1}/{len(factors)}"
        log.info(message)

        # Convert to arrays
        sln = np.array(sln)
        err = np.array(err)
        reg = np.array(reg)

        # Save in a dictionary
        return TikhoTests({"factors": factors, "solution": sln, "error": err, "reg": reg})
