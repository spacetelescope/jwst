"""
Utilities for the ATOCA (Darveau-Bernier 2021, in prep).
ATOCA: Algorithm to Treat Order ContAmination (English)
       Algorithme de Traitement d’Ordres ContAmines (French)

@authors: Antoine Darveau-Bernier, Geert Jan Talens
"""

import numpy as np
from scipy.sparse import find, diags, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d, RectBivariateSpline, Akima1DInterpolator
from scipy.optimize import minimize_scalar, brentq
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# ==============================================================================
# Code for generating indices on the oversampled wavelength grid.
# ==============================================================================


def arange_2d(starts, stops, dtype=None):
    """Create a 2D array containing a series of ranges. The ranges do not have
    to be of equal length.

    Parameters
    ----------
    starts : int or array[int]
        Start values for each range.
    stops : int or array[int]
        End values for each range.
    dtype : str
        Type of the output values.

    Returns
    -------
    out : array[int]
        2D array of ranges.
    mask : array[bool]
        Mask indicating valid elements.
    """

    # Ensure starts and stops are arrays.
    starts = np.asarray(starts)
    stops = np.asarray(stops)

    # Check input for starts and stops is valid.
    if starts.shape != stops.shape and starts.shape != ():
        msg = ('Shapes of starts and stops are not compatible, '
               'they must either have the same shape or starts must be scalar.')
        log.critical(msg)
        raise ValueError(msg)

    if np.any(stops < starts):
        msg = 'stops must be everywhere greater or equal to starts.'
        log.critical(msg)
        raise ValueError(msg)

    # If starts was given as a scalar match its shape to stops.
    if starts.shape == ():
        starts = starts * np.ones_like(stops)

    # Compute the length of each range.
    lengths = (stops - starts).astype(int)

    # Initialize the output arrays.
    nrows = len(stops)
    ncols = np.amax(lengths)
    out = np.ones((nrows, ncols), dtype=dtype)
    mask = np.ones((nrows, ncols), dtype='bool')

    # Compute the indices.
    for irow in range(nrows):
        out[irow, :lengths[irow]] = np.arange(starts[irow], stops[irow])
        mask[irow, :lengths[irow]] = False

    return out, mask


# ==============================================================================
# Code for converting to a sparse matrix and back.
# ==============================================================================


def sparse_k(val, k, n_k):
    """Transform a 2D array `val` to a sparse matrix.

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

    mat = csr_matrix((data, (row, col)), shape=(n_i, n_k))

    return mat


def unsparse(matrix, fill_value=np.nan):
    """Convert a sparse matrix to a 2D array of values and a 2D array of position.

    Parameters
    ----------
    matrix : csr_matrix
        The input sparse matrix.
    fill_value : float
        Value to fill 2D array for undefined positions; default to np.nan

    Returns
    ------
    out : 2d array
        values of the matrix. The shape of the array is given by:
        (matrix.shape[0], maximum number of defined value in a column).
    col_out : 2d array
        position of the columns. Same shape as `out`.
    """

    col, row, val = find(matrix.T)
    n_row, n_col = matrix.shape

    good_rows, counts = np.unique(row, return_counts=True)

    # Define the new position in columns
    i_col = np.indices((n_row, counts.max()))[1]
    i_col = i_col[good_rows]
    i_col = i_col[i_col < counts[:, None]]

    # Create outputs and assign values
    col_out = np.ones((n_row, counts.max()), dtype=int) * -1
    col_out[row, i_col] = col
    out = np.ones((n_row, counts.max())) * fill_value
    out[row, i_col] = val

    return out, col_out


# ==============================================================================
# Code for building wavelength grids.
# ==============================================================================


def get_wave_p_or_m(wave_map, dispersion_axis=1):
    """ Compute upper and lower boundaries of a pixel map,
    given the pixel central value.
    Parameters
    ----------
    wave_map : array[float]
        2d-map of the pixel central wavelength
    dispersion_axis : int, optional
        Which axis is the dispersion axis (0 or 1)

    Returns
    -------
    wave_upper, wave_lower
    The wavelength upper and lower boundaries of each pixel, given the central value.
    """
    # Get wavelength boundaries of each pixels
    wave_left, wave_right = get_wv_map_bounds(wave_map, dispersion_axis=dispersion_axis)

    # The outputs depend on the direction of the spectral axis.
    invalid = (wave_map == 0)
    if ((wave_right >= wave_left) | invalid).all():
        wave_plus, wave_minus = wave_right, wave_left
    elif ((wave_right <= wave_left) | invalid).all():
        wave_plus, wave_minus = wave_left, wave_right
    else:
        msg = 'Some pixels do not follow the expected dispersion axis!'
        log.critical(msg)
        raise ValueError(msg)

    return wave_plus, wave_minus


def get_wv_map_bounds(wave_map, dispersion_axis=1):
    """ Compute boundaries of a pixel map, given the pixel central value.
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
    """
    if dispersion_axis == 1:
        # Simpler to use transpose
        wave_map = wave_map.T
    elif dispersion_axis != 0:
        msg = 'Dispersion axis must be 0 or 1!'
        log.critical(msg)
        raise ValueError(msg)

    # Initialize arrays.
    wave_top = np.zeros_like(wave_map)
    wave_bottom = np.zeros_like(wave_map)

    (n_row, n_col) = wave_map.shape
    for idx in range(n_col):
        wave_col = wave_map[:, idx]

        # Compute the change in wavelength for valid cols
        idx_valid = np.isfinite(wave_col)
        idx_valid &= (wave_col > 0)
        wv_col_valid = wave_col[idx_valid]
        delta_wave = np.diff(wv_col_valid)

        # Init values
        wv_col_top = np.zeros_like(wv_col_valid)
        wv_col_bottom = np.zeros_like(wv_col_valid)

        # Compute the wavelength values on the top and bottom edges of each pixel.
        wv_col_top[1:] = wv_col_valid[:-1] + delta_wave / 2  # TODO check this logic.
        wv_col_top[0] = wv_col_valid[0] - delta_wave[0] / 2
        wv_col_bottom[:-1] = wv_col_valid[:-1] + delta_wave / 2
        wv_col_bottom[-1] = wv_col_valid[-1] + delta_wave[-1] / 2

        wave_top[idx_valid, idx] = wv_col_top
        wave_bottom[idx_valid, idx] = wv_col_bottom

    # De-transpose if it was transposed for computation
    if dispersion_axis == 1:
        wave_top, wave_bottom = wave_top.T, wave_bottom.T

    return wave_top, wave_bottom


def check_dispersion_direction(wave_map, dispersion_axis=1, dwv_sign=-1):
    """Check that the dispersion axis is increasing in the good direction
    given by `dwv_sign``
    Parameters
    ----------
    wave_map : array[float]
        2d-map of the pixel central wavelength
    dispersion_axis : int, optional
        Which axis is the dispersion axis (0 or 1)
    dwv_sign : int, optional
        Direction of increasing wavelengths (-1 or 1)

    Returns
    -------
    bool_map : array[bool]
        Boolean 2d map of the valid dispersion direction, same shape as `wave_map`
    """

    # Estimate the direction of increasing wavelength
    wave_left, wave_right = get_wv_map_bounds(wave_map, dispersion_axis=dispersion_axis)
    dwv = wave_right - wave_left

    # Return bool map of pixels following the good direction
    bool_map = (dwv_sign * dwv >= 0)
    # The bad value could be from left or right so mask both
    bool_map &= np.roll(bool_map, 1, axis=dispersion_axis)

    return bool_map


def mask_bad_dispersion_direction(wave_map, n_max=10, fill_value=0, dispersion_axis=1, dwv_sign=-1):
    """Change value of the pixels in `wave_map` that do not follow the
    general dispersion direction.

    Parameters
    ----------
    wave_map : array[float]
        2d-map of the pixel central wavelength
    n_max : int
        Maximum number of iterations
    fill_value : float
        Value use to replace pixels that do not follow the dispersion direction
    dispersion_axis : int, optional
        Which axis is the dispersion axis (0 or 1)
    dwv_sign : int, optional
        Direction of increasing wavelengths (-1 or 1)

    Returns
    -------
    wave_map : array[float]
        The corrected wave_map.
    convergence flag : bool
        Boolean set to True if all the pixels are now valid, False otherwise.
    """
    # Do not modify the input
    wave_map = wave_map.copy()

    # Make the correction iteratively
    for i_try in range(n_max):
        # Check which pixels are good
        is_good_direction = check_dispersion_direction(wave_map, dispersion_axis, dwv_sign)
        # Stop iteration if all good, or apply correction where needed.
        if is_good_direction.all():
            convergence_flag = True
            break
        else:
            wave_map[~is_good_direction] = fill_value
    else:
        # Did not succeed! :(
        convergence_flag = False

    return wave_map, convergence_flag


def oversample_grid(wave_grid, n_os=1):
    """Create an oversampled version of the input 1D wavelength grid.

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

    # Convert n_os to an array.
    n_os = np.asarray(n_os)

    # n_os needs to have the dimension: len(wave_grid) - 1.
    if n_os.ndim == 0:

        # A scalar was given, repeat the value.
        n_os = np.repeat(n_os, len(wave_grid) - 1)

    elif len(n_os) != (len(wave_grid) - 1):
        # An array of incorrect size was given.
        msg = 'n_os must be a scalar or an array of size len(wave_grid) - 1.'
        log.critical(msg)
        raise ValueError(msg)

    # Grid intervals.
    delta_wave = np.diff(wave_grid)

    # Initialize the new oversampled wavelength grid.
    wave_grid_os = wave_grid.copy()

    # Iterate over oversampling factors to generate new grid points.
    for i_os in range(1, n_os.max()):

        # Consider only intervals that are not complete yet.
        mask = n_os > i_os

        # Compute the new grid points.
        sub_grid = wave_grid[:-1][mask] + (i_os * delta_wave[mask] / n_os[mask])

        # Add the grid points to the oversampled wavelength grid.
        wave_grid_os = np.concatenate([wave_grid_os, sub_grid])

    # Take only unique values and sort them.
    wave_grid_os = np.unique(wave_grid_os)

    return wave_grid_os


def extrapolate_grid(wave_grid, wave_range, poly_ord):
    """Extrapolate the given 1D wavelength grid to cover a given range of values
    by fitting the derivative with a polynomial of a given order and using it to
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

    # Define delta_wave as a function of wavelength by fitting a polynomial.
    delta_wave = np.diff(wave_grid)
    pars = np.polyfit(wave_grid[:-1], delta_wave, poly_ord)
    f_delta = np.poly1d(pars)

    # Extrapolate out-of-bound values on the left-side of the grid.
    grid_left = []
    if wave_range[0] < wave_grid.min():

        # Compute the first extrapolated grid point.
        grid_left = [wave_grid.min() - f_delta(wave_grid.min())]

        # Iterate until the end of wave_range is reached.
        while True:
            next_val = grid_left[-1] - f_delta(grid_left[-1])

            if next_val < wave_range[0]:
                break
            else:
                grid_left.append(next_val)

        # Sort extrapolated vales (and keep only unique).
        grid_left = np.unique(grid_left)

    # Extrapolate out-of-bound values on the right-side of the grid.
    grid_right = []
    if wave_range[-1] > wave_grid.max():

        # Compute the first extrapolated grid point.
        grid_right = [wave_grid.max() + f_delta(wave_grid.max())]

        # Iterate until the end of wave_range is reached.
        while True:
            next_val = grid_right[-1] + f_delta(grid_right[-1])

            if next_val > wave_range[-1]:
                break
            else:
                grid_right.append(next_val)

        # Sort extrapolated vales (and keep only unique).
        grid_right = np.unique(grid_right)

    # Combine the extrapolated sections with the original grid.
    wave_grid_ext = np.concatenate([grid_left, wave_grid, grid_right])

    return wave_grid_ext


def _grid_from_map(wave_map, trace_profile):
    """Define a wavelength grid by taking the wavelength of each column at the
    center of mass of the spatial profile.

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

    # Use only valid columns.
    mask = (trace_profile > 0).any(axis=0) & (wave_map > 0).any(axis=0)

    # Get central wavelength using PSF as weights.
    num = (trace_profile * wave_map).sum(axis=0)
    denom = trace_profile.sum(axis=0)
    center_wv = num[mask] / denom[mask]

    # Make sure the wavelength values are in ascending order.
    sort = np.argsort(center_wv)
    grid = center_wv[sort]

    icols, = np.where(mask)
    return grid, icols[sort]


def grid_from_map(wave_map, trace_profile, wave_range=None, n_os=1, poly_ord=1):
    """Define a wavelength grid by taking the central wavelength at each columns
    given by the center of mass of the spatial profile (so one wavelength per
    column). If wave_range is outside of the wave_map, extrapolate with a
    polynomial of order poly_ord.

    Parameters
    ----------
    wave_map : array[float]
        Array of the pixel wavelengths for a given order.
    trace_profile : array[float]
        Array of the spatial profile for a given order.
    wave_range : list[float]
        Minimum and maximum boundary of the grid to generate, in microns.
        Wave_range must include some wavelengths of wave_map.
    n_os : int or list[int]
        Oversampling of the grid compare to the pixel sampling. Can be
        specified for each order if a list is given. If a single value is given
        it will be used for all orders.
    poly_ord : int
        Order of the polynomial use to extrapolate the grid.

    Returns
    -------
    grid_os : array[float]
        Wavelength grid with oversampling applied
    """

    # Different treatment if wave_range is given.
    if wave_range is None:
        out, _ = _grid_from_map(wave_map, trace_profile)
    else:
        # Get an initial estimate of the grid.
        grid, icols = _grid_from_map(wave_map, trace_profile)

        # Check if extrapolation needed. If so, out_col must be False.
        extrapolate = (wave_range[0] < grid.min()) | (wave_range[1] > grid.max())

        # Make sure grid is between the range
        mask = (wave_range[0] <= grid) & (grid <= wave_range[-1])

        # Check if grid and wv_range are compatible
        if not mask.any():
            msg = "Invalid wave_map or wv_range."
            log.critical(msg)
            raise ValueError(msg)

        grid, icols = grid[mask], icols[mask]

        # Extrapolate values out of the wv_map if needed
        if extrapolate:
            grid = extrapolate_grid(grid, wave_range, poly_ord)

        out = grid

    # Apply oversampling
    grid_os = oversample_grid(out, n_os=n_os)

    return grid_os


def get_soss_grid(wave_maps, trace_profiles, wave_min=0.55, wave_max=3.0, n_os=None):
    """Create a wavelength grid specific to NIRISS SOSS mode observations.
    Assumes 2 orders are given, use grid_from_map if only one order is needed.

    Parameters
    ----------
    wave_maps : array[float]
        Array containing the pixel wavelengths for order 1 and 2.
    trace_profiles : array[float]
        Array containing the spatial profiles for order 1 and 2.
    wave_min : float
        Minimum wavelength the output grid should cover.
    wave_max : float
        Maximum wavelength the output grid should cover.
    n_os : int or list[int]
        Oversampling of the grid compared to the pixel sampling. Can be
        specified for each order if a list is given. If a single value is given
        it will be used for all orders.

    Returns
    -------
    wave_grid_soss : array[float]
        Wavelength grid optimized for extracting SOSS spectra across
        order 1 and order 2.
    """

    # Check n_os input, default value is 2 for all orders.
    if n_os is None:
        n_os = [2, 2]
    elif np.ndim(n_os) == 0:
        n_os = [n_os, n_os]
    elif len(n_os) != 2:
        msg = (f"n_os must be an integer or a 2 element list or array of "
               f"integers, got {n_os} instead")
        log.critical(msg)
        raise ValueError(msg)

    # Generate a wavelength range for each order.
    # Order 1 covers the reddest part of the spectrum,
    # so apply wave_max on order 1 and vice versa for order 2.

    # Take the most restrictive wave_min for order 1
    wave_min_o1 = np.maximum(wave_maps[0].min(), wave_min)

    # Take the most restrictive wave_max for order 2.
    wave_max_o2 = np.minimum(wave_maps[1].max(), wave_max)

    # Now generate range for each orders
    range_list = [[wave_min_o1, wave_max],
                  [wave_min, wave_max_o2]]

    # Use grid_from_map to construct separate oversampled grids for both orders.
    wave_grid_o1 = grid_from_map(wave_maps[0], trace_profiles[0],
                                 wave_range=range_list[0], n_os=n_os[0])
    wave_grid_o2 = grid_from_map(wave_maps[1], trace_profiles[1],
                                 wave_range=range_list[1], n_os=n_os[1])

    # Keep only wavelengths in order 1 that aren't covered by order 2.
    mask = wave_grid_o1 > wave_grid_o2.max()
    wave_grid_o1 = wave_grid_o1[mask]

    # Combine the order 1 and order 2 grids.
    wave_grid_soss = np.concatenate([wave_grid_o1, wave_grid_o2])

    # Sort values (and keep only unique).
    wave_grid_soss = np.unique(wave_grid_soss)

    return wave_grid_soss


def _trim_grids(all_grids, grid_range=None):
    """ Remove all parts of the grids that are not in range
    or that are already covered by grids with higher priority,
    i.e. preceding in the list.
    """
    grids_trimmed = []
    for grid in all_grids:
        # Remove parts of the grid that are not in the wavelength range
        if grid_range is not None:
            # Find where the limit values fall on the grid
            i_min = np.searchsorted(grid, grid_range[0], side='right')
            i_max = np.searchsorted(grid, grid_range[1], side='left')
            # Make sure it is a valid value and take one grid point past the limit
            # since the oversampling could squeeze some nodes near the limits
            i_min = np.max([i_min - 1, 0])
            i_max = np.min([i_max, len(grid) - 1])
            # Trim the grid
            grid = grid[i_min:i_max + 1]

        # Remove parts of the grid that are already covered
        if len(grids_trimmed) > 0:
            # Use all grids already trimmed (so higher in priority)
            conca_grid = np.concatenate(grids_trimmed)
            # Find values below or above
            is_below = grid < np.min(conca_grid)
            is_above = grid > np.max(conca_grid)

            # Do nothing yet if it surrounds the previous grid
            if is_below.any() and is_above.any():
                msg = 'Grid surrounds another grid, better to split in 2 parts.'
                log.warning(msg)

            # Remove values already covered, but keep one
            # index past the limit
            elif is_below.any():
                idx = np.max(np.nonzero(is_below))
                idx = np.min([idx + 1, len(grid) - 1])
                grid = grid[:idx + 1]
            elif is_above.any():
                idx = np.min(np.nonzero(is_above))
                idx = np.max([idx - 1, 0])
                grid = grid[idx:]

            # If all is covered, no need to do it again, so empty grid.
            else:
                grid = np.array([])

        # Save trimmed grid
        grids_trimmed.append(grid)

    return grids_trimmed


def make_combined_adaptive_grid(all_grids, all_estimate, grid_range=None,
                                max_iter=10, rtol=10e-6, tol=0.0, max_total_size=1000000):
    """Return an irregular oversampled grid needed to reach a
    given precision when integrating over each intervals of `grid`.
    The grid is built by subdividing iteratively each intervals that
    did not reach the required precision.
    The precision is computed based on the estimate of the integrals
    using a first order Romberg integration.

    Parameters
    ----------
    all_grid : list[array]
        List of grid (arrays) to pass to adapt_grid, in order of importance.
    all_estimate : list[callable]
        List of function (callable) to estimate the precision needed to oversample the grid.
        Must match the corresponding `grid` in `all_grid`.
    max_iter : int, optional
        Number of times the intervals can be subdivided. The smallest
        subdivison of the grid if max_iter is reached will then be given
        by delta_grid / 2^max_iter. Needs to be greater then zero.
        Default is 10.
    rtol : float, optional
        The desired relative tolerance. Default is 10e-6, so 10 ppm.
    tol : float, optional
        The desired absolute tolerance. Default is 0 to prioritize `rtol`.
    max_total_size : int, optional
        maximum size of the output grid. Default is 1 000 000.
    Returns
    -------
    os_grid : 1D array
        Oversampled combined grid which minimizes the integration error based on
        Romberg's method
    """
    # Save parameters for adapt_grid
    kwargs = dict(max_iter=max_iter, rtol=rtol, tol=tol)

    # Remove unneeded parts of the grids
    all_grids = _trim_grids(all_grids, grid_range=grid_range)

    # Save native size of each grids (use later to adjust max_grid_size)
    all_sizes = [len(grid) for grid in all_grids]

    # Iterate over grids to build the combined grid
    combined_grid = np.array([])  # Init with empty array
    for i_grid, grid in enumerate(all_grids):

        estimate = all_estimate[i_grid]

        # Get the max_grid_size, considering the other grids
        # First, remove length already used
        max_grid_size = max_total_size - combined_grid.size
        # Save some space for next grids (at least the native grid size)
        for i_size, size in enumerate(all_sizes):
            if i_size > i_grid:
                max_grid_size = max_grid_size - size
        # Make sure it is at least the size of the native grid.
        kwargs['max_grid_size'] = np.max([max_grid_size, all_sizes[i_grid]])

        # Oversample the grid based on tolerance required
        grid, is_converged = adapt_grid(grid, estimate, **kwargs)

        # Update grid sizes
        all_sizes[i_grid] = grid.size

        # Check convergence
        if not is_converged:
            msg = 'Precision cannot be garanteed:'
            if grid.size < kwargs['max_grid_size']:
                msg += (f' smallest subdivision 1/{2 ** kwargs["max_iter"]:2.1e}'
                        f' was reached for grid index = {i_grid}')
            else:
                total_size = np.sum(all_sizes)
                msg += ' max grid size of '
                msg += ' + '.join([f'{size}' for size in all_sizes])
                msg += f' = {total_size} was reached for grid index = {i_grid}.'
            log.warning(msg)

        # Remove regions already covered in the output grid
        if len(combined_grid) > 0:
            idx_covered = (np.min(combined_grid) <= grid)
            idx_covered &= (grid <= np.max(combined_grid))
            grid = grid[~idx_covered]

        # Combine grids
        combined_grid = np.concatenate([combined_grid, grid])

    # Sort values (and keep only unique).
    combined_grid = np.unique(combined_grid)

    # Final trim to make sure it respects the range
    if grid_range is not None:
        idx_in_range = (grid_range[0] <= combined_grid)
        idx_in_range &= (combined_grid <= grid_range[-1])
        combined_grid = combined_grid[idx_in_range]

    return combined_grid


def _romberg_diff(b, c, k):
    """Compute the differences for the Romberg quadrature corrections.
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

    tmp = 4.0**k
    diff = (tmp * c - b) / (tmp - 1.0)

    return diff


def _difftrap(fct, intervals, numtraps):
    """Perform part of the trapezoidal rule to integrate a function. Assume that
    we had called difftrap with all lower powers-of-2 starting with 1. Calling
    difftrap only returns the summation of the new ordinates. It does not
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


def get_n_nodes(grid, fct, divmax=10, tol=1.48e-4, rtol=1.48e-4):
    """Refine parts of a grid to reach a specified integration precision
    based on Romberg integration of a callable function or method.
    Returns the number of nodes needed in each intervals of
    the input grid to reach the specified tolerance over the integral
    of `fct` (a function of one variable).

    Note: This function is based on scipy.integrate.quadrature.romberg. The
    difference between it and the scipy version is that it is vectorized to deal
    with multiple intervals separately. It also returns the number of nodes
    needed to reached the required precision instead of returning the value of
    the integral.

    Parameters
    ----------
    grid : array[float]
        Grid for integration. Each section of this grid is treated as a
        separate integral; if grid has length N, N-1 integrals are optimized.
    fct : callable
        Function to be integrated.
    divmax : int
        Maximum order of extrapolation.
    tol : float
        The desired absolute tolerance.
    rtol : float
        The desired relative tolerance.

    Returns
    -------
    n_grid : array[int]
        Number of nodes needed on each distinct intervals in the grid to reach
        the specified tolerance.
    residual : array[float]
        Estimate of the error in each intervals. Same length as n_grid.
    """

    # Initialize some variables.
    n_intervals = len(grid) - 1
    i_bad = np.arange(n_intervals)
    n_grid = np.repeat(-1, n_intervals)
    residual = np.repeat(np.nan, n_intervals)

    # Change the 1D grid into a 2D set of intervals.
    intervals = np.array([grid[:-1], grid[1:]])
    intrange = np.diff(grid)
    err = np.inf

    # First estimate without subdivision.
    numtraps = 1
    ordsum = _difftrap(fct, intervals, numtraps)
    results = intrange * ordsum
    last_row = [results]

    for i_div in range(1, divmax + 1):

        # Increase the number of trapezoids by factors of 2.
        numtraps *= 2

        # Evaluate trapz integration for intervals that are not converged.
        ordsum += _difftrap(fct, intervals[:, i_bad], numtraps)
        row = [intrange[i_bad] * ordsum / numtraps]

        # Compute Romberg for each of the computed sub grids.
        for k in range(i_div):
            romb_k = _romberg_diff(last_row[k], row[k], k + 1)
            row = np.vstack([row, romb_k])

        # Save R(n,n) and R(n-1, n-1) from Romberg method.
        results = row[i_div]
        lastresults = last_row[i_div - 1]

        # Estimate error.
        err = np.abs(results - lastresults)

        # Find intervals that are converged.
        conv = (err < tol) | (err < rtol * np.abs(results))

        # Save number of nodes for these intervals.
        n_grid[i_bad[conv]] = numtraps

        # Save residuals.
        residual[i_bad] = err

        # Stop if all intervals have converged.
        if conv.all():
            break

        # Find intervals not converged.
        i_bad = i_bad[~conv]

        # Save last_row and ordsum for the next iteration for non-converged
        # intervals.
        ordsum = ordsum[~conv]
        last_row = row[:, ~conv]

    else:
        # Warn that convergence is not reached everywhere.
        log.warning(f"divmax {divmax} exceeded. Latest difference = {err.max()}")

    # Make sure all values of n_grid where assigned during the process.
    if (n_grid == -1).any():
        msg = f"Values where not assigned at grid position: {np.where(n_grid == -1)}"
        log.critical(msg)
        raise ValueError(msg)

    return n_grid, residual


def estim_integration_err(grid, fct):
    """Estimate the integration error on each intervals
    of the grid using 1rst order Romberg integration.

    Parameters
    ----------
    grid: 1d array [float]
        Grid for integration. Each sections of this grid are treated
        as separate integrals. So if grid has length N; N-1 integrals are
        tested.
    fct: callable
        Function to be integrated.

    Returns
    -------
    err, rel_err: error and relative error of each integrations, with length = length(grid) - 1
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
    non_zero = (romb != 0)
    rel_err = np.full_like(err, np.inf)
    rel_err[non_zero] = np.abs(err[non_zero] / romb[non_zero])

    return err, rel_err


def adapt_grid(grid, fct, max_iter=10, rtol=10e-6, tol=0.0, max_grid_size=None):
    """Return an irregular oversampled grid needed to reach a
    given precision when integrating over each intervals of `grid`.
    The grid is built by subdividing iteratively each intervals that
    did not reach the required precision.
    The precision is computed based on the estimate of the integrals
    using a first order Romberg integration.

    Parameters
    ----------
    grid: array
        Grid for integration. Each sections of this grid are treated
        as separate integrals. So if grid has length N; N-1 integrals are
        optimized.
    fct: callable
        Function to be integrated. Must be a function of `grid`
    max_iter: int, optional
        Number of times the intervals can be subdivided. The smallest
        subdivison of the grid if max_iter is reached will then be given
        by delta_grid / 2^max_iter. Needs to be greater then zero.
        Default is 10.
    rtol: float, optional
        The desired relative tolerance. Default is 10e-6, so 10 ppm.
    tol: float, optional
        The desired absolute tolerance. Default is 0 to prioritize `rtol`.
    max_grid_size: int, optional
        maximum size of the output grid. Default is None, so no constraint.
    Returns
    -------
    os_grid  : 1D array
        Oversampled grid which minimizes the integration error based on
        Romberg's method
    convergence_flag: bool
        Whether the estimated tolerance was reach everywhere or not.
    See Also
    --------
    scipy.integrate.quadrature.romberg
    References
    ----------
    [1] 'Romberg's method' https://en.wikipedia.org/wiki/Romberg%27s_method

    """
    # No limit of max_grid_size not given
    if max_grid_size is None:
        max_grid_size = np.inf

    # Init some flags
    max_size_reached = (grid.size >= max_grid_size)

    # Iterate until precision is reached of max_iter
    for _ in range(max_iter):

        # Estimate error using Romberg integration
        err, rel_err = estim_integration_err(grid, fct)

        # Check where precision is reached
        converged = (err < tol) | (rel_err < rtol)
        is_converged = converged.all()

        # Check if max grid size was reached
        if max_size_reached or is_converged:
            # Then stop iteration
            break

        # Intervals that didn't reach the precision will be subdivided
        n_oversample = np.full(err.shape, 2, dtype=int)
        # No subdivision for the converged ones
        n_oversample[converged] = 1

        # Check if the maximum size will be reached.
        # If so, prioritize the intervals with the largest estimated errors
        # to reach the maximum size
        os_grid_size = n_oversample.sum()
        if os_grid_size > max_grid_size:
            # How many nodes can be added to reach max?
            n_nodes_remaining = max_grid_size - grid.size

            # Find the position of the nodes with the largest error
            idx_largest_err = np.argsort(rel_err)[-n_nodes_remaining:]

            # Build new oversample array and assign only largest errors
            n_oversample = np.ones(err.shape, dtype=int)
            n_oversample[idx_largest_err] = 2

            # Flag to stop iterations
            max_size_reached = True

        # Generate oversampled grid (subdivide)
        grid = oversample_grid(grid, n_os=n_oversample)

        # Make sure sorted and unique.
        grid = np.unique(grid)

    return grid, is_converged


# ==============================================================================
# Code for handling the throughput and kernels.
# ==============================================================================


class ThroughputSOSS(interp1d):

    def __init__(self, wavelength, throughput):
        """Create an instance of scipy.interpolate.interp1d to handle the
        throughput values.

        Parameters
        ----------
        wavelength : array[float]
            A wavelength array.
        throughput : array[float]
            The throughput values corresponding to the wavelengths.
        """

        # Interpolate
        super().__init__(wavelength, throughput, kind='cubic', fill_value=0,
                         bounds_error=False)


class WebbKernel:  # TODO could probably be cleaned-up somewhat, may need further adjustment.

    def __init__(self, wave_kernels, kernels, wave_map, n_os, n_pix,  # TODO kernels may need to be flipped?
                 bounds_error=False, fill_value="extrapolate"):
        """A handler for the kernel values.

        Parameters
        ----------
        wave_kernels : array[float]
            Kernels for wavelength array.
        kernels : array[float]
            Kernels for throughput array.
        wave_map : array[float]
            Wavelength map of the detector. Since WebbPSF returns kernels in
            the pixel space, we need a wave_map to convert to wavelength space.
        n_os : int
            Oversampling of the kernels.
        n_pix : int
            Length of the kernels in pixels.
        bounds_error : bool
            If True, raise an error when trying to call the function out of the
            interpolation range. If False, the values will be extrapolated.
        fill_value : str
            How to extrapolate when needed. Only default "extrapolate"
            currently implemented.
        """

        # Mask where wv_map is equal to 0
        wave_map = np.ma.array(wave_map, mask=(wave_map == 0))

        # Force wv_map to have the red wavelengths
        # at the end of the detector
        if np.diff(wave_map, axis=-1).mean() < 0:
            wave_map = np.flip(wave_map, axis=-1)

        # Number of columns
        ncols = wave_map.shape[-1]

        # Create oversampled pixel position array
        pixels = np.arange(-(n_pix // 2), n_pix // 2 + (1 / n_os), (1 / n_os))

        # `wave_kernel` has only the value of the central wavelength
        # of the kernel at each points because it's a function
        # of the pixels (so depends on wv solution).
        wave_center = wave_kernels[0, :]

        # Use the wavelength solution to create a mapping between pixels and wavelengths
        # First find the all kernels that fall on the detector.
        wave_min = np.amin(wave_map[wave_map > 0])
        wave_max = np.amax(wave_map[wave_map > 0])
        i_min = np.searchsorted(wave_center, wave_min)
        i_max = np.searchsorted(wave_center, wave_max) - 1

        # Use the next kernels at each extremities to define the
        # boundaries of the interpolation to use in the class
        # RectBivariateSpline (at the end)
        bbox = [None, None,
                wave_center[np.maximum(i_min - 1, 0)],
                wave_center[np.minimum(i_max + 1, len(wave_center) - 1)]]
        #######################

        # Keep only kernels that fall on the detector.
        kernels = kernels[:, i_min:i_max + 1].copy()
        wave_kernels = wave_kernels[:, i_min:i_max + 1].copy()
        wave_center = np.array(wave_kernels[0, :])

        # Save minimum kernel value (greater than zero)
        kernels_min = np.min(kernels[(kernels > 0.0)])

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
            i_row, i_col = np.unravel_index(
                np.argmin(np.abs(wave_map - wv_c)), wave_map.shape
            )
            # Update wavelength center value
            # (take the nearest pixel center value)
            wave_center[i_cen] = wave_map[i_row, i_col]

            # Surrounding columns
            index = i_col + i_surround

            # Make sure it's on the detector
            i_good = (index >= 0) & (index < ncols)

            # Assign wv values
            wv[i_good] = wave_map[i_row, index[i_good]]

            # Fit n=1 polynomial
            poly_i = np.polyfit(i_surround[~wv.mask], wv[~wv.mask], 1)

            # Project on os pixel grid
            wave_kernels[:, i_cen] = np.poly1d(poly_i)(pixels)

            # Save coeffs
            poly.append(poly_i)

        # Save attributes
        self.n_pix = n_pix
        self.n_os = n_os
        self.wave_kernels = wave_kernels
        self.kernels = kernels
        self.pixels = pixels
        self.wave_center = wave_center
        self.poly = np.array(poly)
        self.fill_value = fill_value
        self.bounds_error = bounds_error
        self.min_value = kernels_min

        # 2d Interpolate
        self.f_ker = RectBivariateSpline(pixels, wave_center, kernels, bbox=bbox)

    def __call__(self, wave, wave_c):
        """Returns the kernel value, given the wavelength and the kernel central
         wavelength.

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
        fill_value = self.fill_value
        bounds_error = self.bounds_error
        n_wv_c = len(wave_center)
        f_ker = self.f_ker
        n_pix = self.n_pix
        min_value = self.min_value

        # #################################
        # First, convert wv value in pixels
        # using a linear interpolation
        # #################################

        # Find corresponding interval
        i_wv_c = np.searchsorted(wave_center, wave_c) - 1

        # Deal with values out of bounds
        if bounds_error:
            message = "Value of wv center out of interpolation range"
            log.critical(message)
            raise ValueError(message)
        elif fill_value == "extrapolate":
            i_wv_c[i_wv_c < 0] = 0
            i_wv_c[i_wv_c >= (n_wv_c - 1)] = n_wv_c - 2
        else:
            message = f"`fill_value`={fill_value} is not an valid option."
            log.critical(message)
            raise ValueError(message)

        # Compute coefficients that interpolate along wv_centers
        d_wv_c = wave_center[i_wv_c + 1] - wave_center[i_wv_c]
        a_c = (wave_center[i_wv_c + 1] - wave_c) / d_wv_c
        b_c = (wave_c - wave_center[i_wv_c]) / d_wv_c

        # Compute a_pix and b_pix from the equation:
        # pix = a_pix * lambda + b_pix
        a_pix = 1 / (a_c * poly[i_wv_c, 0] + b_c * poly[i_wv_c + 1, 0])
        b_pix = -(a_c * poly[i_wv_c, 1] + b_c * poly[i_wv_c + 1, 1])
        b_pix /= (a_c * poly[i_wv_c, 0] + b_c * poly[i_wv_c + 1, 0])

        # Compute pixel values
        pix = a_pix * wave + b_pix

        # ######################################
        # Second, compute kernel value on the
        # interpolation grid (pixel x wv_center)
        # ######################################

        webbker = f_ker(pix, wave_c, grid=False)

        # Make sure it's not negative and greater than the min value
        webbker = np.clip(webbker, min_value, None)

        # and put out-of-range values to zero.
        webbker[pix > n_pix // 2] = 0
        webbker[pix < -(n_pix // 2)] = 0

        return webbker


# ==============================================================================
# Code for building the convolution matrix (c matrix).
# ==============================================================================


def gaussians(x, x0, sig, amp=None):
    """Gaussian function

    Parameters
    ----------
    x : array[float]
        Array of points over which gaussian to be defined.
    x0 : float
        Center of the gaussian.
    sig : float
        Standard deviation of the gaussian.
    amp : float
        Value of the gaussian at the center.

    Returns
    -------
    values : array[float]
        Array of gaussian values for input x.
    """

    # Amplitude term
    if amp is None:
        amp = 1. / np.sqrt(2. * np.pi * sig**2.)

    values = amp * np.exp(-0.5 * ((x - x0) / sig) ** 2.)

    return values


def fwhm2sigma(fwhm):
    """Convert a full width half max to a standard deviation, assuming a gaussian

    Parameters
    ----------
    fwhm : float
        Full-width half-max of a gaussian.

    Returns
    -------
    sigma : float
        Standard deviation of a gaussian.
    """

    sigma = fwhm / np.sqrt(8. * np.log(2.))

    return sigma


def to_2d(kernel, grid_range):
    """ Build a 2d kernel array with a constant 1D kernel (input)

    Parameters
    ----------
    kernel : array[float]
        Input 1D kernel.
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
    kernel_2d = np.tile(kernel, (n_k_c, 1)).T

    return kernel_2d


def _get_wings(fct, grid, h_len, i_a, i_b):
    """Compute values of the kernel at grid[+-h_len]

    Parameters
    ----------
    fct : callable
        Function that returns the value of the kernel, given
        a grid value and the center of the kernel.
        fct(grid, center) = kernel
        grid and center have the same length.
    grid : array[float]
        grid where the kernel is projected
    h_len : int
        Half-length where we compute kernel value.
    i_a : int
        Index of grid axis 0 where to apply convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[i_a:i_b].
    i_b : int
        index of grid axis 1 where to apply convolution.

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
    grid_new = grid[i_grid:i_b - h_len]

    # Re-use dummy variable `i_grid`
    i_grid = len(grid_new)

    # Compute kernel at the left end.
    # `i_grid` accounts for smaller length.
    ker = fct(grid_new, grid[i_b - i_grid:i_b])
    left[-i_grid:] = ker

    # Add the right value on the grid
    # Possibility that it falls out of the grid;
    # take last value of the grid if so.
    # Same steps as the left end (see above)
    i_grid = np.min([n_k, i_b + h_len])
    grid_new = grid[i_a + h_len:i_grid]
    i_grid = len(grid_new)
    ker = fct(grid_new, grid[i_a:i_a + i_grid])
    right[:i_grid] = ker

    return left, right


def trpz_weight(grid, length, shape, i_a, i_b):
    """Compute weights due to trapezoidal integration

    Parameters
    ----------
    grid : array[float]
        grid where the integration is projected
    length : int
        length of the kernel
    shape : tuple[int]
        shape of the compact convolution 2d array
    i_a : int
        Index of grid axis 0 where to apply convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[i_a:i_b].
    i_b : int
        index of grid axis 1 where to apply convolution.

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


def fct_to_array(fct, grid, grid_range, thresh=1e-5, length=None):
    """Build a compact kernel 2d array based on a kernel function
    and a grid to project the kernel

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
    thresh : float, optional
        Threshold to cut the kernel wings. If `length` is specified,
        `thresh` will be ignored.
    length : int, optional
        Length of the kernel. Must be odd.

    Returns
    -------
    kern_array : array[float]
        2D array of kernel projected onto grid.
    """

    # Assign range where the convolution is defined on the grid
    i_a, i_b = grid_range

    # Init with the value at kernel's center
    out = fct(grid, grid)[i_a:i_b]

    # Add wings
    if length is None:
        # Generate a 2D array of the grid iteratively until
        # thresh is reached everywhere.

        # Init parameters
        length = 1
        h_len = 0  # Half length

        # Add value on each sides until thresh is reached
        while True:
            # Already update half-length
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
                left[left < thresh] = 0.
                right[right < thresh] = 0.

                # add new values to output
                out = np.vstack([left, out, right])

        # Weights due to integration (from the convolution)
        weights = trpz_weight(grid, length, out.shape, i_a, i_b)

    elif (length % 2) == 1:  # length needs to be odd
        # Generate a 2D array of the grid iteratively until
        # specified length is reached.

        # Compute number of half-length
        n_h_len = (length - 1) // 2

        # Simply iterate to compute needed wings
        for h_len in range(1, n_h_len + 1):
            # Compute next left and right ends of the kernel
            left, right = _get_wings(fct, grid, h_len, i_a, i_b)

            # Add new kernel values
            out = np.vstack([left, out, right])

        # Weights due to integration (from the convolution)
        weights = trpz_weight(grid, length, out.shape, i_a, i_b)

    else:
        msg = "`length` provided to `fct_to_array` must be odd."
        log.critical(msg)
        raise ValueError(msg)

    kern_array = (out * weights)
    return kern_array


def cut_ker(ker, n_out=None, thresh=None):
    """Apply a cut on the convolution matrix boundaries.

    Parameters
    ----------
    ker : array[float]
        convolution kernel in compact form, so
        shape = (N_ker, N_k_convolved)
    n_out : int, list[int] or tuple[int]
        Number of kernel's grid point to keep on the boundaries.
        If an int is given, the same number of points will be
        kept on each boundaries of the kernel (left and right).
        If 2 elements are given, it corresponds to the left and right
        boundaries.
    thresh : float
        threshold used to determine the boundaries cut.
        If n_out is specified, this is ignored.

    Returns
    ------
    ker : array[float]
        The same kernel matrix as the input ker, but with the cut applied.
    """

    # Assign kernel length and number of kernels
    n_ker, n_k_c = ker.shape

    # Assign half-length of the kernel
    h_len = (n_ker - 1) // 2

    # Determine n_out with thresh if not given
    if n_out is None:

        if thresh is None:
            # No cut to apply
            return ker
        else:
            # Find where to cut the kernel according to thresh
            i_left = np.where(ker[:, 0] >= thresh)[0][0]
            i_right = np.where(ker[:, -1] >= thresh)[0][-1]

            # Make sure it is on the good wing. Take center if not.
            i_left = np.minimum(i_left, h_len)
            i_right = np.maximum(i_right, h_len)

    # Else, unpack n_out
    else:
        # Could be a scalar or a 2-elements object)
        try:
            i_left, i_right = n_out
        except TypeError:
            i_left, i_right = n_out, n_out

        # Find the position where to cut the kernel
        # Make sure it is not out of the kernel grid,
        # so i_left >= 0 and i_right <= len(kernel)
        i_left = np.maximum(h_len - i_left, 0)
        i_right = np.minimum(h_len + i_right, n_ker - 1)

    # Apply the cut
    for i_k in range(0, i_left):
        # Add condition in case the kernel is larger
        # than the grid where it's projected.
        if i_k < n_k_c:
            ker[:i_left - i_k, i_k] = 0

    for i_k in range(i_right + 1 - n_ker, 0):
        # Add condition in case the kernel is larger
        # than the grid where it's projected.
        if -i_k <= n_k_c:
            ker[i_right - n_ker - i_k:, i_k] = 0

    return ker


def sparse_c(ker, n_k, i_zero=0):
    """Convert a convolution kernel in compact form (N_ker, N_k_convolved)
    to sparse form (N_k_convolved, N_k)

    Parameters
    ----------
    ker : array[float]
        Convolution kernel in compact form, with shape (N_kernel, N_kc)
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
        err_msg = "Length of the convolution kernel given to sparse_c should be odd."
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
    matrix = diags(diag_val, offset, shape=(n_k_c, n_k), format="csr")

    return matrix


def get_c_matrix(kernel, grid, bounds=None, i_bounds=None, norm=True,
                 sparse=True, n_out=None, thresh_out=None, **kwargs):
    """Return a convolution matrix
    Can return a sparse matrix (N_k_convolved, N_k)
    or a matrix in the compact form (N_ker, N_k_convolved).
    N_k is the length of the grid on which the convolution
    will be applied, N_k_convolved is the length of the
    grid after convolution and N_ker is the maximum length of
    the kernel. If the default sparse matrix option is chosen,
    the convolution can be applied on an array f | f = fct(grid)
    by a simple matrix multiplication:
    f_convolved = c_matrix.dot(f)

    Parameters
    ----------
    kernel: ndarray (1D or 2D), callable
        Convolution kernel. Can be already 2D (N_ker, N_k_convolved),
        giving the kernel for each items of the convolved grid.
        Can be 1D (N_ker), so the kernel is the same. Can be a callable
        with the form f(x, x0) where x0 is the position of the center of
        the kernel. Must return a 1D array (len(x)), so a kernel value
        for each pairs of (x, x0). If kernel is callable, the additional
        kwargs `thresh` and `length` will be used to project the kernel.
    grid: one-d-array:
        The grid on which the convolution will be applied.
        For example, if C is the convolution matrix,
        f_convolved = C.f(grid)
    bounds: 2-elements object
        The bounds of the grid on which the convolution is defined.
        For example, if bounds = (a,b),
        then grid_convolved = grid[a <= grid <= b].
        It dictates also the dimension of f_convolved
    sparse: bool, optional
        return a sparse matrix (N_k_convolved, N_k) if True.
        return a matrix (N_ker, N_k_convolved) if False.
    n_out: integer or 2-integer object, optional
        Specify how to deal with the ends of the convolved grid.
        `n_out` points will be used outside from the convolved
        grid. Can be different for each ends if 2-elements are given.
    thresh_out: float, optional
        Specify how to deal with the ends of the convolved grid.
        Points with a kernel value less then `thresh_out` will
        not be used outside from the convolved grid.
    thresh: float, optional
        Only used when `kernel` is callable to define the maximum
        length of the kernel. Truncate when `kernel` < `thresh`
    length: int, optional
        Only used when `kernel` is callable to define the maximum
        length of the kernel.
    """

    # Define range where the convolution is defined on the grid.
    # If `i_bounds` is not specified, try with `bounds`.
    if i_bounds is None:

        if bounds is None:
            a, b = 0, len(grid)
        else:
            a = np.min(np.where(grid >= bounds[0])[0])
            b = np.max(np.where(grid <= bounds[1])[0]) + 1

    else:
        # Make sure it is absolute index, not relative
        # So no negative index.
        if i_bounds[1] < 0:
            i_bounds[1] = len(grid) + i_bounds[1]

        a, b = i_bounds

    # Generate a 2D kernel depending on the input
    if callable(kernel):
        kernel = fct_to_array(kernel, grid, [a, b], **kwargs)
    elif kernel.ndim == 1:
        kernel = to_2d(kernel, [a, b])

    if kernel.ndim != 2:
        msg = ("Input kernel to get_c_matrix must be callable or"
               " array with one or two dimensions.")
        log.critical(msg)
        raise ValueError(msg)
    # Kernel should now be a 2-D array (N_kernel x N_kc)

    # Normalize if specified
    if norm:
        kernel = kernel / np.nansum(kernel, axis=0)

    # Apply cut for kernel at boundaries
    kernel = cut_ker(kernel, n_out, thresh_out)

    if sparse:
        # Convert to a sparse matrix.
        kernel = sparse_c(kernel, len(grid), a)

    return kernel


class NyquistKer:
    """Define a gaussian convolution kernel at the nyquist
    sampling. For a given point on the grid x_i, the kernel
    is given by a gaussian with
    FWHM = n_sampling * (dx_(i-1) + dx_i) / 2.
    The FWHM is computed for each elements of the grid except
    the extremities (not defined). We can then generate FWHM as
    a function of the grid and interpolate/extrapolate to get
    the kernel as a function of its position relative to the grid.
    """

    def __init__(self, grid, n_sampling=2, bounds_error=False,
                 fill_value="extrapolate", **kwargs):
        """Parameters
        ----------
        grid : array[float]
            Grid used to define the kernels
        n_sampling : int, optional
            Sampling of the grid.
        bounds_error : bool
            Argument for `interp1d` to get FWHM as a function of the grid.
        fill_value : str
            Argument for `interp1d` to choose fill method to get FWHM.
        """

        # Delta grid
        d_grid = np.diff(grid)

        # The full width half max is n_sampling
        # times the mean of d_grid
        fwhm = (d_grid[:-1] + d_grid[1:]) / 2
        fwhm *= n_sampling

        # What we really want is sigma, not FWHM
        sig = fwhm2sigma(fwhm)

        # Now put sigma as a function of the grid
        sig = interp1d(grid[1:-1], sig, bounds_error=bounds_error,
                       fill_value=fill_value, **kwargs)

        self.fct_sig = sig

    def __call__(self, x, x0):
        """Parameters
        ----------
        x : array[float]
            position where the kernel is evaluated
        x0 : array[float]
            position of the kernel center for each x.

        Returns
        -------
        Value of the gaussian kernel for each set of (x, x0)
        """

        # Get the sigma of each gaussian
        sig = self.fct_sig(x0)

        return gaussians(x, x0, sig)


# ==============================================================================
# Code for doing Tikhonov regularisation.
# ==============================================================================


def finite_diff(x):
    """Returns the finite difference matrix operator based on x.

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

    # Build matrix
    diff_matrix = diags([-1.], shape=(n_x - 1, n_x))
    diff_matrix += diags([1.], 1, shape=(n_x - 1, n_x))

    return diff_matrix


def finite_second_d(grid):
    """Returns the second derivative operator based on grid

    Parameters
    ----------
    grid : array[float]
        1D array where the second derivative will be computed.

    Returns
    -------
    second_d : array[float]
        Operator to compute the second derivative, so that
        f" = second_d.dot(f), where f is a function
        projected on `grid`.
    """

    # Finite difference operator
    d_matrix = finite_diff(grid)

    # Delta lambda
    d_grid = d_matrix.dot(grid)

    # First derivative operator
    first_d = diags(1. / d_grid).dot(d_matrix)

    # Second derivative operator
    second_d = finite_diff(grid[:-1]).dot(first_d)

    # don't forget the delta lambda
    second_d = diags(1. / d_grid[:-1]).dot(second_d)

    return second_d


def finite_first_d(grid):
    """Returns the first derivative operator based on grid

    Parameters
    ----------
    grid : array[float]
        Grid where the first derivative will be computed.

    Returns
    -------
    first_d : array[float]
        Operator to compute the second derivative, so that
        f' = first_d.dot(f), where f is a function
        projected on `grid`.
    """

    # Finite difference operator
    d_matrix = finite_diff(grid)

    # Delta lambda
    d_grid = d_matrix.dot(grid)

    # First derivative operator
    first_d = diags(1. / d_grid).dot(d_matrix)

    return first_d


def get_tikho_matrix(grid, n_derivative=1, d_grid=True, estimate=None, pwr_law=0):
    """Wrapper to return the tikhonov matrix given a grid and the derivative degree.

    Parameters
    ----------
    grid : array[float]
        1D grid where the Tikhonov matrix is projected
    n_derivative : int, optional
        Degree of derivative. Possible values are 1 or 2
    d_grid : bool, optional
        Whether to divide the differential operator by the grid differences,
        which corresponds to an actual approximation of the derivative or not.
    estimate : callable (preferably scipy.interpolate.UnivariateSpline), optional
        Estimate of the solution on which the tikhonov matrix is applied.
        Must be a function of `grid`. If UnivariateSpline, then the derivatives
        are given directly (so best option), otherwise the tikhonov matrix will be
        applied to `estimate(grid)`. Note that it is better to use `d_grid=True`
    pwr_law: float, optional
        Power law applied to the scale differentiated estimate, so the estimate
        of tikhonov_matrix.dot(solution). It will be applied as follows:
        norm_factor * scale_factor.dot(tikhonov_matrix)
        where scale_factor = 1/(estimate_derivative)**pwr_law
        and norm_factor = 1/sum(scale_factor)
    Returns
    -------
    t_mat : array[float]
        The tikhonov matrix.
    """
    if d_grid:
        input_grid = grid
    else:
        input_grid = np.arange(len(grid))

    if n_derivative == 1:
        t_mat = finite_first_d(input_grid)
    elif n_derivative == 2:
        t_mat = finite_second_d(input_grid)
    else:
        msg = "`n_derivative` must be 1 or 2."
        log.critical(msg)
        raise ValueError(msg)

    if estimate is not None:
        if hasattr(estimate, 'derivative'):
            # Get the derivatives directly from the spline
            if n_derivative == 1:
                derivative = estimate.derivative(n=n_derivative)
                tikho_factor_scale = derivative(grid[:-1])
            elif n_derivative == 2:
                derivative = estimate.derivative(n=n_derivative)
                tikho_factor_scale = derivative(grid[1:-1])
        else:
            # Apply tikho matrix on estimate
            tikho_factor_scale = t_mat.dot(estimate(grid))

        # Make sure all positive
        tikho_factor_scale = np.abs(tikho_factor_scale)
        # Apply power law
        # (similar to 'kunasz1973'?)
        tikho_factor_scale = np.power(tikho_factor_scale, -pwr_law)
        # Normalize
        valid = np.isfinite(tikho_factor_scale)
        tikho_factor_scale /= np.sum(tikho_factor_scale[valid])

        # If some values are not finite, set to the max value
        # so it will be more regularized
        valid = np.isfinite(tikho_factor_scale)
        if not valid.all():
            value = np.max(tikho_factor_scale[valid])
            tikho_factor_scale[~valid] = value

        # Apply to tikhonov matrix
        t_mat = diags(tikho_factor_scale).dot(t_mat)

    return t_mat


def curvature_finite(factors, log_reg2, log_chi2):
    """Compute the curvature in log space using finite differences

    Parameters
    ----------
    factors : array[float]
        Regularisation factors (not in log).
    log_reg2 : array[float]
        norm-2 of the regularisation term (in log10).
    log_chi2 : array[float]
        norm-2 of the chi2 term (in log10).

    Returns
    -------
    factors : array[float]
        Sorted and cut version of input factors array.
    curvature : array[float]

    """
    # Make sure it is sorted according to the factors
    idx = np.argsort(factors)
    factors, log_chi2, log_reg2 = factors[idx], log_chi2[idx], log_reg2[idx]

    # Get first and second derivatives
    chi2_deriv = get_finite_derivatives(factors, log_chi2)
    reg2_deriv = get_finite_derivatives(factors, log_reg2)

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


def get_finite_derivatives(x_array, y_array):
    """ Compute first and second finite derivatives
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
    """ Generate array given the relative range around an index.

    Parameters
    ----------
    idx : int
        Center index value
    relative_range : iterable[int]
        relative bounds around center value to create new array
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
    out = np.arange(*abs_range, 1)

    return out


def _minimize_on_grid(factors, val_to_minimize, interpolate, interp_index=None):
    """ Find minimum of a grid using akima spline interpolation to get a finer estimate

    Parameters
    ----------
    factors : array[float]
        1D array of Tikhonov factors for which value array is calculated
    val_to_minimize : array[float]
        1D array of values to be minimized, e.g. chi^2 or curvature.
    interpolate: bool, optional
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
        opt_args = {"bounds": bounds,
                    "method": "bounded"}
        min_fac = minimize_scalar(fct, **opt_args).x

        # Back to linear scale
        min_fac = 10. ** min_fac

    else:
        # Simply return the min value
        # if no interpolation required
        min_fac = factors[idx_min]

    return min_fac


def _find_intersect(factors, y_val, thresh, interpolate, search_range=None):
    """ Find the root of y_val - thresh (so the intersection between thresh and y_val)
    Parameters
    ----------
    factors : array[float]
        1D array of Tikhonov factors for which value array is calculated
    y_val : array[float]
        1D array of values.
    thresh: float
        Threshold use in 'd_chi2' mode. Find the highest factor where the
        derivative of the chi2 derivative is below thresh.
    interpolate: bool, optional
        If True, use interpolation to find a finer minimum;
        otherwise, return minimum value in array.
    search_range : iterable[int], optional
        Relative range of grid indices around the value to interpolate.
        If not specified, defaults to [0,3].

    Returns
    -------
    float
        Factor corresponding to the best approximation of the intersection
        point.
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
    cond_below = (y_val < thresh)
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
        d_chi2_spl = interp1d(x_val, y_val - thresh, kind='linear')

        # Use index only around the best value
        max_length = len(y_val)
        index = _get_interp_idx_array(idx_below, search_range, max_length)

        # Find the root
        bracket = (x_val[index[0]], x_val[index[-1]])
        best_val = brentq(d_chi2_spl, *bracket)

        # Back to linear scale
        best_val = 10. ** best_val

    else:
        # Simply return the value
        best_val = factors[idx_below]

    return best_val


def soft_l1(z):
    return 2 * ((1 + z)**0.5 - 1)


def cauchy(z):
    return np.log(1 + z)


def linear(z):
    return z


LOSS_FUNCTIONS = {'soft_l1': soft_l1, 'cauchy': cauchy, 'linear': linear}


class TikhoTests(dict):
    """
    Class to save Tikhonov tests for different factors.
    All the tests are stored in the attribute `tests` as a dictionary

    Parameters
    ----------
    test_dict : dict
        Dictionary holding arrays for `factors`, `solution`, `error`, and `reg`
        by default.
    """

    DEFAULT_TRESH_DERIVATIVE = (('chi2', 1e-5),
                                ('chi2_soft_l1', 1e-4),
                                ('chi2_cauchy', 1e-3))

    def __init__(self, test_dict=None, default_chi2='chi2_cauchy'):
        """
        Parameters
        ----------
        test_dict : dict
            Dictionary holding arrays for `factors`, `solution`, `error`, and `reg`
            by default.
        default_chi2: string
            Type of chi2 loss used by default. Options are chi2, chi2_soft_l1, chi2_cauchy.
        """
        # Define the number of data points
        # (length of the "b" vector in the tikhonov regularisation)
        if test_dict is None:
            print('Unable to get the number of data points. Setting `n_points` to 1')
            n_points = 1
        else:
            n_points = len(test_dict['error'][0].squeeze())

        # Save attributes
        self.n_points = n_points
        self.default_chi2 = default_chi2
        self.default_thresh = {chi2_type: thresh
                               for (chi2_type, thresh)
                               in self.DEFAULT_TRESH_DERIVATIVE}

        # Initialize so it behaves like a dictionary
        super().__init__(test_dict)

        chi2_loss = {'chi2': 'linear',
                     'chi2_soft_l1': 'soft_l1',
                     'chi2_cauchy': 'cauchy'}
        for chi2_type, loss in chi2_loss.items():
            try:
                # Save the chi2
                self[chi2_type]
            except KeyError:
                self[chi2_type] = self.compute_chi2(loss=loss)
#         # Save different loss function for chi2
#         self['chi2_soft_l1'] = self.compute_chi2(loss='soft_l1')
#         self['chi2_cauchy'] = self.compute_chi2(loss='cauchy')

    def compute_chi2(self, tests=None, n_points=None, loss='linear'):
        """ Calculates the reduced chi squared statistic

        Parameters
        ----------
        tests : dict, optional
            Dictionary from which we take the error array; if not provided,
            self is used
        n_points : int, optional
            Number of data points; if not provided, self.n_points is used

        Returns
        -------
        float
            Sum of the squared error array divided by the number of data points
        """
        # If not given, take the tests from the object
        if tests is None:
            tests = self

        # Get the loss function
        if isinstance(loss, str):
            try:
                loss = LOSS_FUNCTIONS[loss]
            except KeyError as e:
                keys = [key for key in LOSS_FUNCTIONS.keys()]
                msg = f'loss={loss} not a valid key. Must be one of {keys} or callable.'
                raise e(msg)
        elif not callable(loss):
            raise ValueError('Invalid value for loss.')

        # Compute the reduced chi^2 for all tests
        chi2 = np.nanmean(loss(tests['error']**2), axis=-1)
        # Remove residual dimensions
        chi2 = chi2.squeeze()

        return chi2

    def get_chi2_derivative(self, key=None):
        """ Compute derivative of the chi2 with respect to log10(factors)

        Parameters
        ----------
        key: str
            which chi2 is used for computations. Default is self.default_chi2.

        Returns
        -------
        factors_leftd : array[float]
            factors array, shortened to match length of derivative.
        d_chi2 : array[float]
            derivative of chi squared array with respect to log10(factors)
        """
        if key is None:
            key = self.default_chi2

        # Compute finite derivative
        fac_log = np.log10(self['factors'])
        d_chi2 = np.diff(self[key]) / np.diff(fac_log)

        # Update size of factors to fit derivatives
        # Equivalent to derivative on the left side of the nodes
        factors_leftd = self['factors'][1:]

        return factors_leftd, d_chi2

    def compute_curvature(self, tests=None, key=None):

        if key is None:
            key = self.default_chi2

        # If not given, take the tests from the object
        if tests is None:
            tests = self

        # Compute the curvature...
        # Get the norm-2 of the regularisation term
        reg2 = np.nansum(tests['reg'] ** 2, axis=-1)

        factors, curv = curvature_finite(tests['factors'],
                                         np.log10(self[key]),
                                         np.log10(reg2))

        return factors, curv

    def best_tikho_factor(self, tests=None, interpolate=True, interp_index=None,
                          mode='curvature', key=None, thresh=None):
        """Compute the best scale factor for Tikhonov regularisation.
        It is determined by taking the factor giving the highest logL on
        the detector or the highest curvature of the l-curve,
        depending on the chosen mode.
        Parameters
        ----------
        tests : dictionary, optional
            Results of tikhonov extraction tests for different factors.
            Must have the keys "factors" and "-logl". If not specified,
            the tests from self.tikho.tests are used.
        interpolate : bool, optional
            If True, use spline interpolation to find a finer minimum.
            Default is true.
        interp_index : list, optional
            Relative range of grid indices around the minimum value to
            interpolate across. If not specified, defaults to [-2,4].
        mode : string
            How to find the best factor: 'chi2', 'curvature' or 'd_chi2'.
        thresh : float
            Threshold for use in 'd_chi2' mode. Find the highest factor where
            the derivative of the chi2 derivative is below thresh.

        Returns
        -------
        float
            Best scale factor as determined by the selected algorithm
        """
        if key is None:
            key = self.default_chi2

        if thresh is None:
            thresh = self.default_thresh[key]

        # Use pre-run tests if not specified
        if tests is None:
            tests = self

        # Number of factors
        n_fac = len(tests['factors'])

        # Determine the mode (what do we minimize?)
        if mode == 'curvature' and n_fac > 2:
            # Compute the curvature
            factors, curv = tests.compute_curvature()

            # Find min factor
            best_fac = _minimize_on_grid(factors, curv, interpolate, interp_index)

        elif mode == 'chi2':
            # Simply take the chi2 and factors
            factors = tests['factors']
            y_val = tests[key]

            # Find min factor
            best_fac = _minimize_on_grid(factors, y_val, interpolate, interp_index)

        elif mode == 'd_chi2' and n_fac > 1:
            # Compute the derivative of the chi2
            factors, y_val = tests.get_chi2_derivative()

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
            best_fac = _find_intersect(factors[idx], y_val[idx], thresh, interpolate, interp_index)

        elif mode in ['curvature', 'd_chi2', 'chi2']:
            best_fac = np.max(tests['factors'])
            msg = (f'Could not compute {mode} because number of factor={n_fac}. '
                   f'Setting best factor to max factor: {best_fac:.5e}')
            log.warning(msg)

        else:
            msg = (f'`mode`={mode} is not a valid option for '
                   f'TikhoTests.best_tikho_factor().')
            log.critical(msg)
            raise ValueError(msg)

        # Return estimated best scale factor
        return best_fac


class Tikhonov:
    """
    Tikhonov regularization to solve the ill-posed problem A.x = b, where
    A is accidentally singular or close to singularity. Tikhonov regularization
    adds a regularization term in the equation and aim to minimize the
    equation: ||A.x - b||^2 + ||gamma.x||^2
    where gamma is the Tikhonov regularization matrix.
    """

    def __init__(self, a_mat, b_vec, t_mat, valid=True):
        """
        Parameters
        ----------
        a_mat : matrix-like object (2d)
            matrix A in the system to solve A.x = b
        b_vec : vector-like object (1d)
            vector b in the system to solve A.x = b
        t_mat : matrix-like object (2d)
            Tikhonov regularisation matrix to be applied on b_vec.
        valid : bool, optional
            If True, solve the system only for valid indices. The
            invalid values will be set to np.nan. Default is True.
        """

        # Save input matrix
        self.a_mat = a_mat
        self.b_vec = b_vec
        self.t_mat = t_mat

        # Pre-compute some matrix for the linear system to solve
        t_mat_2 = (t_mat.T).dot(t_mat)  # squared tikhonov matrix
        a_mat_2 = a_mat.T.dot(a_mat)  # squared model matrix
        result = (a_mat.T).dot(b_vec.T)
        idx_valid = (result.toarray() != 0).squeeze()  # valid indices to use if `valid` is True

        # Save pre-computed matrix
        self.t_mat_2 = t_mat_2
        self.a_mat_2 = a_mat_2
        self.result = result
        self.idx_valid = idx_valid

        # Save other attributes
        self.valid = valid
        self.test = None

        return

    def solve(self, factor=1.0):
        """
        Minimize the equation ||A.x - b||^2 + ||gamma.x||^2
        by solving (A_T.A + gamma_T.gamma).x = A_T.b
        gamma is the Tikhonov matrix multiplied by a scale factor

        Parameters
        ----------
        factor : float, optional
            multiplicative constant of the regularization matrix

        Returns
        ------
        array
            Solution of the system (1d array)
        """
        # Get needed attributes
        a_mat_2 = self.a_mat_2
        result = self.result
        t_mat_2 = self.t_mat_2
        valid = self.valid
        idx_valid = self.idx_valid

        # Matrix gamma squared (with scale factor)
        gamma_2 = factor ** 2 * t_mat_2

        # Finalize building matrix
        matrix = a_mat_2 + gamma_2

        # Initialize solution
        solution = np.full(matrix.shape[0], np.nan)

        # Consider only valid indices if in valid mode
        if valid:
            idx = idx_valid
        else:
            idx = np.full(len(solution), True)

        # Solve
        matrix = matrix[idx, :][:, idx]
        result = result[idx]
        solution[idx] = spsolve(matrix, result)

        return solution

    def test_factors(self, factors):
        """
        Test multiple factors

        Parameters
        ----------
        factors : array[float]
            1D array of factors to test

        Returns
        ------
        dict
            Dictionary of test results
        """

        log.info('Testing factors...')

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
            err.append(a_mat.dot(sln[-1]) - b_vec)

            # Save regularization term
            reg_i = t_mat.dot(sln[-1])
            reg.append(reg_i)

            # Print
            message = '{}/{}'.format(i_fac, len(factors))
            log.info(message)

        # Final print
        message = '{}/{}'.format(i_fac + 1, len(factors))
        log.info(message)

        # Convert to arrays
        sln = np.array(sln)
        err = np.array(err)
        reg = np.array(reg)

        # Save in a dictionary

        tests = TikhoTests({'factors': factors,
                            'solution': sln,
                            'error': err,
                            'reg': reg})

        return tests
