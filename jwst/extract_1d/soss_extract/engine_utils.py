#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# General imports.
import numpy as np
from warnings import warn
from scipy.integrate import AccuracyWarning
from scipy.sparse import find, diags, identity, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d, RectBivariateSpline

# Plotting.
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# ==============================================================================
# Code for generating indices on the oversampled wavelength grid.
# ==============================================================================


def arange_2d(starts, stops, dtype=None):
    """Create a 2D array containing a series of ranges. The ranges do not have
    to be of equal length.

    :param starts: start values for each range.
    :param stops: end values for each range.
    :param dtype: the type of the output values.

    :type starts: int or array[int]
    :type stops: int or array[int]
    :type dtype: str

    :returns: out, mask - 2D array of ranges and a mask indicating valid
        elements.
    :rtype: Tuple(array[int], array[bool])
    """

    # Ensure starts and stops are arrays.
    starts = np.asarray(starts)
    stops = np.asarray(stops)

    # Check input for starts and stops is valid.
    if (starts.shape != stops.shape) & (starts.shape != ()):
        msg = ('Shapes of starts and stops are not compatible, '
               'they must either have the same shape or starts must be scalar.')
        raise ValueError(msg)

    if np.any(stops < starts):
        msg = 'stops must be everywhere greater or equal to starts.'
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
    """
    Transform a 2D array `val` to a sparse matrix.
    `k` is use for the position in the second axis
    of the matrix. The resulting sparse matrix will
    have the shape : ((len(k), n_k))
    Set k elements to a negative value when not defined
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
    """
    Convert a sparse matrix to a 2D array of values and a 2D array of position.

    Returns
    ------
    out: 2d array
        values of the matrix. The shape of the array is given by:
        (matrix.shape[0], maximum number of defined value in a column).
    col_out: 2d array
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


def get_wave_p_or_m(wave_map):
    # TODO rename function?
    """Compute lambda_plus and lambda_minus of pixel map, given the pixel
    central value.

    :param wave_map: Array of the pixel wavelengths for a given order.

    :type wave_map: array[float]

    :returns: wave_plus, wave_minus - The wavelength edges of each pixel,
        given the central value.
    :rtype: Tuple(array[float], array[float])
    """

    wave_map = wave_map.T  # Simpler to use transpose

    # Iniitialize arrays.
    wave_left = np.zeros_like(wave_map)
    wave_right = np.zeros_like(wave_map)

    # Compute the change in wavelength.
    delta_wave = np.diff(wave_map, axis=0)

    # Compute the wavelength values on the left and right edges of each pixel.
    wave_left[1:] = wave_map[:-1] + delta_wave/2  # TODO check this logic.
    wave_left[0] = wave_map[0] - delta_wave[0]/2
    wave_right[:-1] = wave_map[:-1] + delta_wave/2
    wave_right[-1] = wave_map[-1] + delta_wave[-1]/2

    # The outputs depend on the direction of the spectral axis.
    if (wave_right >= wave_left).all():
        wave_plus, wave_minus = wave_right.T, wave_left.T
    elif (wave_right <= wave_left).all():
        wave_plus, wave_minus = wave_left.T, wave_right.T
    else:
        raise ValueError('Bad pixel values for wavelength.')

    return wave_plus, wave_minus


def oversample_grid(wave_grid, n_os=1):
    """Create an oversampled version of the input 1D wavelength grid.

    :param wave_grid: Wavelength grid to be oversampled.
    :param n_os: Oversampling factor. If it is a scalar, take the same value for each
        interval of the grid. If it is an array, n_os specifies the oversampling
        at each interval of the grid, so len(n_os) = len(wave_grid) - 1.

    :type wave_grid: array[float]
    :type n_os: int or array[int]

    :returns: wave_grid_os - The oversampled wavelength grid.
    :rtype: array[float]
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
        sub_grid = (wave_grid[:-1][mask] + i_os*delta_wave[mask]/n_os[mask])

        # Add the grid points to the oversampled wavelength grid.
        wave_grid_os = np.concatenate([wave_grid_os, sub_grid])

    # Take only uniqyue values and sort them.
    wave_grid_os = np.unique(wave_grid_os)

    return wave_grid_os


def extrapolate_grid(wave_grid, wave_range, poly_ord):
    """Extrapolate the given 1D wavelength grid to cover a given range of values
    by fitting the derivate with a polynomial of a given order and using it to
    compute subsequent values at both ends of the grid.

    :param wave_grid: Wavelength grid to be extrapolated.
    :param wave_range: Wavelength range the new grid should cover.
    :param poly_ord: Order of the polynomial used to fit the derivative of
        wave_grid.

    :type wave_grid: array[float]
    :type wave_range: list[float]
    :type poly_ord: int

    :returns: wave_grid_ext - The extrapolated 1D wavelength grid.
    :rtype: array[float]
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


def _grid_from_map(wave_map, aperture, out_col=False):
    # TODO is out_col still needed.
    """Define a wavelength grid by taking the wavelength of each column at the
    center of mass of the spatial profile.

    :param wave_map: Array of the pixel wavelengths for a given order.
    :param aperture: Array of the spatial profile for a given order.
    :param out_col:

    :type wave_map: array[float]
    :type aperture: array[float]
    :type out_col: bool

    :returns:
    :rtype:
    """

    # Use only valid columns.
    mask = (aperture > 0).any(axis=0) & (wave_map > 0).any(axis=0)

    # Get central wavelength using PSF as weights.
    num = (aperture * wave_map).sum(axis=0)
    denom = aperture.sum(axis=0)
    center_wv = num[mask]/denom[mask]

    # Make sure the wavelength values are in ascending order.
    sort = np.argsort(center_wv)
    grid = center_wv[sort]

    if out_col:  # TODO I don't like this type of contruction much.
        # Return index of columns if out_col is True.
        icols, = np.where(mask)
        return grid, icols[sort]
    else:
        # Return sorted and unique if out_col is False.
        grid = np.unique(grid)
        return grid


def grid_from_map(wave_map, aperture, wave_range=None, n_os=1, poly_ord=1,
                  out_col=False):
    # TODO is out_col still needed.
    """Define a wavelength grid by taking the central wavelength at each columns
    given by the center of mass of the spatial profile (so one wavelength per
    column). If wave_range is outside of the wave_map, extrapolate with a
    polynomial of order poly_ord.

    :param wave_map: Array of the pixel wavelengths for a given order.
    :param aperture: Array of the spatial profile for a given order.
    :param wave_range: Minimum and maximum boundary of the grid to generate,
        in microns. wave_range must include some wavelenghts of wave_map.
    :param n_os: Oversampling of the grid compare to the pixel sampling. Can be
        specified for each order if a list is given. If a single value is given
        it will be used for all orders.
    :param poly_ord: Order of the polynomial use to extrapolate the grid.
        Default is 1.
    :param out_col: Return columns. TODO It will be forced to False if extrapolation is needed.

    :type wave_map: array[float]
    :type aperture: array[float]
    :type wave_range: List[float]
    :type poly_ord: int
    :type out_col: bool
    :type n_os: int

    :returns:
    :rtype:
    """

    # Different treatement if wave_range is given.
    if wave_range is None:
        out = _grid_from_map(wave_map, aperture, out_col=out_col)
    else:
        # Get an initial estimate of the grid.
        grid, icols = _grid_from_map(wave_map, aperture, out_col=True)

        # Check if extrapolation needed. If so, out_col must be False.
        extrapolate = (wave_range[0] < grid.min()) | (wave_range[1] > grid.max())
        if extrapolate and out_col:
            out_col = False
            msg = ("Cannot extrapolate and return columns. "
                   "Setting out_col = False.")
            warn(msg)

        # Make sure grid is between the range
        mask = (wave_range[0] <= grid) & (grid <= wave_range[-1])

        # Check if grid and wv_range are compatible
        if not mask.any():
            msg = "Invalid wave_map or wv_range. wv_range: {}"
            raise ValueError(msg.format(wave_range))

        grid, icols = grid[mask], icols[mask]

        # Extrapolate values out of the wv_map if needed
        if extrapolate:
            grid = extrapolate_grid(grid, wave_range, poly_ord)

        # Different output depending on `out_col`
        if out_col:
            out = grid, icols
        else:
            out = grid

    # Apply oversampling
    if out_col:
        # Return grid and columns TODO this doesn't seem to do that? Would crash if out was a tuple?
        return [oversample_grid(out_i, n_os=n_os) for out_i in out]
    else:
        # Only the grid
        return oversample_grid(out, n_os=n_os)


def get_soss_grid(wave_maps, apertures, wave_min=0.55, wave_max=3.0, n_os=None):
    """Create a wavelength grid specific to NIRISS SOSS mode observations.
    Assumes 2 orders are given, use grid_from_map if only one order is needed.

    Parameters
    ----------
    :param wave_maps: Array containing the pixel wavelengths for order 1 and 2.
    :param apertures: Array containing the spatial profiles for order 1 and 2.
    :param wave_min: Minimum wavelength the output grid should cover.
    :param wave_max: Maximum wavelength the output grid should cover.
    :param n_os: Oversampling of the grid compare to the pixel sampling. Can be
        specified for each order if a list is given. If a single value is given
        it will be used for all orders.

    :type wave_maps: array[float]
    :type apertures: array[float]
    :type wave_min: float
    :type wave_max: float
    :type n_os: int or List[int]

    :returns: wave_grid_soss - A wavelength grid optimized for extracting SOSS
        spectra across order 1 and order 2.
    :rtype: array[float]
    """

    # Check n_os input, default value is 2 for all orders.
    if n_os is None:
        n_os = [2, 2]
    elif np.ndim(n_os) == 0:
        n_os = [n_os, n_os]
    elif len(n_os) != 2:
        msg = ("n_os must be an integer or a 2 element list or array of "
               "integers, got {} instead")
        raise ValueError(msg.format(n_os))

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
    wave_grid_o1 = grid_from_map(wave_maps[0], apertures[0],
                                 wave_range=range_list[0], n_os=n_os[0])
    wave_grid_o2 = grid_from_map(wave_maps[1], apertures[1],
                                 wave_range=range_list[1], n_os=n_os[1])

    # Keep only wavelengths in order 1 that aren't covered by order 2.
    mask = wave_grid_o1 > wave_grid_o2.max()
    wave_grid_o1 = wave_grid_o1[mask]

    # Combine the order 1 and order 2 grids.
    wave_grid_soss = np.concatenate([wave_grid_o1, wave_grid_o2])

    # Sort values (and keep only unique).
    wave_grid_soss = np.unique(wave_grid_soss)

    return wave_grid_soss


def _romberg_diff(b, c, k):
    """Compute the differences for the Romberg quadrature corrections.
    See Forman Acton's "Real Computing Made Real," p 143.

    :param b: R(n-1, m-1) of Rombergs method.
    :param c: R(n, m-1) of Rombergs method.
    :param k: The parameter m of Rombergs method.

    :type b: float or array[float]
    :type c: float or array[float]
    :type k: int

    :returns: R(n, m) of Rombergs method.
    :rtype: float or array[float]
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

    :param fct: Function to be integrated.
    :param intervals: A 2D array of integration intervals of shape (Nx2) or a
        single interval of shape (2,).
    :param numtraps: The number of trapezoids used to integrate the interval.
        numtraps must be a power of 2.

    :type fct: callable
    :type intervals: array[float]
    :type numtraps: int

    :returns: s - The sum of function values at the new trapezoid boundaries
    compared to numtraps = numtraps/2. When numtraps = 1 they are fivided by
    two.
    :rtype: float
    """

    # Convert input intervals to numpy array
    intervals = np.asarray(intervals)

    # If intervals is 1D assume it's a single interval.
    if intervals.ndim == 1:
        intervals = intervals[:, np.newaxis]

    # Check the value of numtraps.
    if numtraps <= 0:  # TODO check it is a power of 2? Or change input to log2(numtraps)?
        raise ValueError("numtraps must be > 0 in difftrap().")

    if numtraps == 1:
        # Return the function evaluations for a single trapezoid.
        # Only points add the edge of the interval need to be halfed.
        ordsum = 0.5*(fct(intervals[0]) + fct(intervals[1]))

    else:
        # Number of new points compared to lower 2**N multiple of trapezoids.
        numtosum = numtraps/2

        # Find coordinates of new points.
        h = (intervals[1] - intervals[0])/numtosum
        lox = intervals[0] + 0.5*h
        points = lox[np.newaxis, :] + h*np.arange(numtosum)[:, np.newaxis]

        # Evalaute and sum the new points.
        ordsum = np.sum(fct(points), axis=0)

    return ordsum


def get_n_nodes(grid, fct, divmax=10, tol=1.48e-4, rtol=1.48e-4):
    """Refine parts of a grid to reach a specified integration precision
    based on Romberg integration of a callable function or method.
    Returns the number of nodes needed in each intervals of
    the input grid to reach the specified tolerance over the integral
    of `fct` (a function of one variable).

    Note: This function is based on scipy.integrate.quadrature.romberg. The
    difference between it and the scipy version is that it is vectorised to deal
    with multiple intervals separately. It also returns the number of nodes
    needed to reached the required precision instead of returning the value of
    the integral.

    :param grid: Grid for integration. Each sections of this grid are treated
        as separate integrals. So if grid has length N; N-1 integrals are
        optimized.
    :param fct: Function to be integrated.
    :param tol: The desired absolute tolerance. Default is 1.48e-4.
    :param rtol: The desired relative tolerance. Default is 1.48e-4.
    :param divmax: Maximum order of extrapolation. Default is 10.

    :type grid: array[float]
    :type fct: callable
    :type tol: float
    :type rtol: float
    :type divmax: int

    :returns: n_grid - Number of nodes needed on each distinct intervals in the
        grid to reach the specified tolerance. If out_res=True also returns
        residual - Estimate of the error in each intervals. Same length as
        n_grid.
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

        # Evaluate trpz integration for intervals that are not converged.
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
        msg = "divmax {%d} exceeded. Latest difference = {}"
        warn(msg.format(divmax, err.max()), AccuracyWarning)

    # Make sure all values of n_grid where assigned during the process.
    if (n_grid == -1).any():
        msg = "Values where not assigned at grid position: {}"
        raise ValueError(msg.format(np.where(n_grid == -1)))

    return n_grid, residual


# ==============================================================================
# Code for handling the throughput and kernels.
# ==============================================================================


class ThroughputSOSS(interp1d):

    def __init__(self, wavelength, throughput):
        """Create an instance of scipy.interpolate.interp1d to handle the
        throughput values.

        :param wavelength: A wavelength array.
        :param throughput: The throughput values corresponding to the
            wavelengths.

        :type wavelength: array[float]
        :type throughput: array[float]
        """

        # Interpolate
        super().__init__(wavelength, throughput, kind='cubic', fill_value=0,
                         bounds_error=False)


class WebbKernel:  # TODO could probably be cleaned-up somewhat, may need further adjustment.

    def __init__(self, wave_kernels, kernels, wave_map, n_os, n_pix,  # TODO kernels may need to be flipped?
                 bounds_error=False, fill_value="extrapolate"):
        """A handler for the kernel values.

        :param wave_kernels:
        :param kernels:
        :param wave_map: Wavelength map of the detector. Since WebbPSF returns
            kernels in the pixel space, we need a wave_map to convert to
            wavelength space.
        :param n_os: Oversampling of the kernels.
        :param n_pix: Length of the kernels in pixels.
        :param bounds_error: If True, raise an error when trying to call the
            function out of the interpolation range. If False, the values will
            be extrapolated. Default is False
        :param fill_value: How to extrapolate when needed. Default is
            "extrapolate". No other options have been implemented.

        :type wave_kernels: array[float]
        :type kernels: array[float]
        :type wave_map: array[float]
        :type n_os: int
        :type n_pix: int
        :type bounds_error: bool
        :type fill_value: str
        """

        # Mask where wv_map is equal to 0
        wave_map = np.ma.array(wave_map, mask=(wave_map == 0))

        # Force wv_map to have the red wavelengths
        # at the end of the detector
        if np.diff(wave_map, axis=-1).mean() < 0:
            wave_map = np.flip(wave_map, axis=-1)

        # Number of columns
        ncols = wave_map.shape[-1]

        # Create oversampled pixel position array  # TODO easier to read form?
        pixels = np.arange(-(n_pix//2), n_pix//2 + 1/n_os, 1/n_os)

        # `wave_kernel` has only the value of the central wavelength
        # of the kernel at each points because it's a function
        # of the pixels (so depends on wv solution).
        wave_center = wave_kernels[0, :]

        # Use the wavelength solution to create a mapping.
        # First find the kernels that fall on the detector.
        wave_min = np.amin(wave_map[wave_map > 0])
        wave_max = np.amax(wave_map[wave_map > 0])
        i_min = np.searchsorted(wave_center, wave_min)  # TODO searchsorted has offsets?
        i_max = np.searchsorted(wave_center, wave_max) - 1

        # SAVE FOR LATER ###########
        # Use the next kernels at each extremities to define the
        # boundaries of the interpolation to use in the class
        # RectBivariateSpline (at the end)
        # bbox = [min pixel, max pixel, min wv_center, max wv_center]
        bbox = [None, None,
                wave_center[np.maximum(i_min-1, 0)],
                wave_center[np.minimum(i_max+1, len(wave_center)-1)]]
        #######################

        # Keep only kernels that fall on the detector.
        kernels, wave_kernels = kernels[:, i_min:i_max+1], wave_kernels[:, i_min:i_max+1]
        wave_center = np.array(wave_kernels[0, :])

        # Then find the pixel closest to each kernel center
        # and use the surrounding pixels (columns)
        # to get the wavelength. At the boundaries,
        # wavelenght might not be defined or falls out of
        # the detector, so fit a 1-order polynomial to
        # extrapolate. The polynomial is also used to interpolate
        # for oversampling.
        i_surround = np.arange(-(n_pix//2), n_pix//2 + 1)
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

        # 2d Interpolate
        self.f_ker = RectBivariateSpline(pixels, wave_center, kernels, bbox=bbox)

    def __call__(self, wave, wave_c):
        """Returns the kernel value, given the wavelength and the kernel central
         wavelength.

        :param wave: Wavelength where the kernel is projected.
        :param wave_c: Central wavelength of the kernel.

        :type wave: array[float]
        :type wave_c: array[float]

        :returns: out - The kernel value.
        """

        wave_center = self.wave_center
        poly = self.poly
        fill_value = self.fill_value
        bounds_error = self.bounds_error
        n_wv_c = len(wave_center)
        f_ker = self.f_ker
        n_pix = self.n_pix

        # #################################
        # First, convert wv value in pixels
        # using a linear interpolation
        # #################################

        # Find corresponding interval
        i_wv_c = np.searchsorted(wave_center, wave_c) - 1

        # Deal with values out of bounds
        if bounds_error:
            message = "Value of wv center out of interpolation range"
            raise ValueError(message)
        elif fill_value == "extrapolate":
            i_wv_c[i_wv_c < 0] = 0
            i_wv_c[i_wv_c >= (n_wv_c - 1)] = n_wv_c - 2
        else:
            message = "`fill_value`={} is not an valid option."
            raise ValueError(message.format(fill_value))

        # Compute coefficients that interpolate along wv_centers
        d_wv_c = wave_center[i_wv_c + 1] - wave_center[i_wv_c]
        a_c = (wave_center[i_wv_c + 1] - wave_c) / d_wv_c
        b_c = (wave_c - wave_center[i_wv_c]) / d_wv_c

        # Compute a_pix and b_pix from the equation:
        # pix = a_pix * lambda + b_pix
        a_pix = 1 / (a_c * poly[i_wv_c, 0] + b_c * poly[i_wv_c+1, 0])
        b_pix = -(a_c * poly[i_wv_c, 1] + b_c * poly[i_wv_c+1, 1])
        b_pix /= (a_c * poly[i_wv_c, 0] + b_c * poly[i_wv_c+1, 0])

        # Compute pixel values
        pix = a_pix * wave + b_pix

        # ######################################
        # Second, compute kernel value on the
        # interpolation grid (pixel x wv_center)
        # ######################################

        webbker = f_ker(pix, wave_c, grid=False)

        # Make sure it's not negative, and put out of range values to zero.
        webbker[webbker < 0] = 0
        webbker[pix > n_pix//2] = 0
        webbker[pix < -(n_pix//2)] = 0

        return webbker

    def show(self):
        """
        Plot kernels.
        The first figure is a 2d image of the kernels.
        The second figure is a 1d image of the kernels
        in the wavelength space.
        """

        # 2D figure of the kernels
        fig1 = plt.figure()

        # Log plot, so clip values <= 0
        image = np.clip(self.kernels, np.min(self.kernels[self.kernels > 0]), np.inf)

        # plot
        plt.pcolormesh(self.wave_center, self.pixels,  image, norm=LogNorm())

        # Labels and others
        plt.colorbar(label="Kernel")
        plt.ylabel("Position relative to center [pixel]")
        plt.xlabel(r"Center wavelength [$\mu m$]")
        plt.tight_layout()

        # 1D figure of all kernels
        fig2 = plt.figure()
        plt.plot(self.wave_kernels, self.kernels)

        # Labels and others
        plt.ylabel("Kernel")
        plt.xlabel(r"Wavelength [$\mu m$]")
        plt.tight_layout()

        return fig1, fig2


# ==============================================================================
# Code for building the convolution matrix (c matrix).
# ==============================================================================


def gaussians(x, x0, sig, amp=None):
    """
    Gaussian function
    """

    # Amplitude term
    if amp is None:
        amp = 1/np.sqrt(2 * np.pi * sig**2)

    values = amp * np.exp(-0.5 * ((x - x0) / sig) ** 2)

    return values


def fwhm2sigma(fwhm):
    """
    Convert a full width half max to a standard deviation, assuming a gaussian
    """

    sigma = fwhm / np.sqrt(8 * np.log(2))

    return sigma


def to_2d(kernel, grid, grid_range):  # TODO parameter grid unused, better name kernel_2d?
    """ Build a 2d kernel array with a constant 1D kernel (input) """

    # Assign range where the convolution is defined on the grid
    a, b = grid_range

    # Get length of the convolved axis
    n_k_c = b - a

    # Return a 2D array with this length
    kernel_2d = np.tile(kernel, (n_k_c, 1)).T

    return kernel_2d


def _get_wings(fct, grid, h_len, i_a, i_b):
    """
    Compute values of the kernel at grid[+-h_len]

    Parameters
    ----------
    fct: callable
        Function that returns the value of the kernel, given
        a grid value and the center of the kernel.
        fct(grid, center) = kernel
        grid and center have the same length.
    grid: 1d array
        grid where the kernel is projected
    i_a, i_b: int
        index of the grid where to apply the convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[i_a:i_b].
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
    i_grid = np.max([0, i_a-h_len])

    # Save the new grid
    grid_new = grid[i_grid:i_b-h_len]

    # Re-use dummy variable `i_grid`
    i_grid = len(grid_new)

    # Compute kernel at the left end.
    # `i_grid` accounts for smaller length.
    ker = fct(grid_new, grid[i_b-i_grid:i_b])
    left[-i_grid:] = ker

    # Add the right value on the grid
    # Possibility that it falls out of the grid;
    # take last value of the grid if so.
    # Same steps as the left end (see above)
    i_grid = np.min([n_k, i_b + h_len])
    grid_new = grid[i_a+h_len:i_grid]
    i_grid = len(grid_new)
    ker = fct(grid_new, grid[i_a:i_a+i_grid])
    right[:i_grid] = ker

    return left, right


def trpz_weight(grid, length, shape, i_a, i_b):
    """
    Compute weights due to trapeze integration

    Parameters
    ----------
    grid: 1d array
        grid where the integration is proojected
    length: int
        length of the kernel
    shape: 2-elements tuple
        shape of the compact convolution 2d array
    i_a, i_b: int
        index of the grid where to apply the convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[i_a:i_b].
    """

    # Index of each element on the convolution matrix
    # with respect to the non-convolved grid
    # `i_grid` has the shape (N_k_convolved, kernel_length - 1)
    i_grid = np.indices(shape)[0] - (length // 2)
    i_grid = np.arange(i_a, i_b)[None, :] + i_grid[:-1, :]

    # Set values out of grid to -1
    i_bad = (i_grid < 0) | (i_grid >= len(grid)-1)
    i_grid[i_bad] = -1

    # Delta lambda
    d_grid = np.diff(grid)

    # Compute weights from trapezoidal integration
    weight = 1/2 * d_grid[i_grid]
    weight[i_bad] = 0

    # Fill output
    out = np.zeros(shape)
    out[:-1] += weight
    out[1:] += weight

    return out


def fct_to_array(fct, grid, grid_range, thresh=1e-5, length=None):
    """
    Build a compact kernel 2d array based on a kernel function
    and a grid to project the kernel

    Parameters
    ----------
    fct: callable
        Function that returns the value of the kernel, given
        a grid value and the center of the kernel.
        fct(grid, center) = kernel
        grid and center have the same length.
    grid: 1d array
        grid where the kernel is projected
    grid_range: 2 element list or tuple
        index of the grid where to apply the convolution.
        Once the convolution applied, the convolved grid will be
        equal to grid[grid_range[0]:grid_range[1]].
    thresh: float, optional
        threshold to cut the kernel wings. If `length` is specified,
        `thresh` will not be taken into account.
    length: int, optional
        length of the kernel. Needs to be odd.
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
        for h_len in range(1, n_h_len+1):
            # Compute next left and right ends of the kernel
            left, right = _get_wings(fct, grid, h_len, i_a, i_b)

            # Add new kernel values
            out = np.vstack([left, out, right])

        # Weights due to integration (from the convolution)
        weights = trpz_weight(grid, length, out.shape, i_a, i_b)

    else:
        raise ValueError("`length` must be odd.")

    return out*weights


def cut_ker(ker, n_out=None, thresh=None):
    """
    Apply a cut on the convolution matrix boundaries.

    Parameters
    ----------
    ker: 2d array
        convolution kernel in compact form, so
        shape = (N_ker, N_k_convolved)
    n_out: int or 2-element int object (list, tuple, etc.)
        Number of kernel's grid point to keep on the boundaries.
        If an int is given, the same number of points will be
        kept on each boundaries of the kernel (left and right).
        If 2 elements are given, it corresponds to the left and right
        boundaries.
    thresh: float
        threshold used to determine the boundaries cut.
        If n_out is specified, this has no effect.

    Returns
    ------
    the same kernel matrix has the input ker, but with the cut applied.
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
            ker[:i_left-i_k, i_k] = 0

    for i_k in range(i_right + 1 - n_ker, 0):
        # Add condition in case the kernel is larger
        # than the grid where it's projected.
        if -i_k <= n_k_c:
            ker[i_right-n_ker-i_k:, i_k] = 0

    return ker


def sparse_c(ker, n_k, i_zero=0):
    """
    Convert a convolution kernel in compact form (N_ker, N_k_convolved)
    to sparse form (N_k_convolved, N_k)

    Parameters
    ----------
    ker : 2d array, (N_kernel, N_kc)
        Convolution kernel in compact form.
    n_k: int
        length of the original grid
    i_zero: int
        position of the first element of the convolved grid
        in the original grid.
    """

    # Assign kernel length and convolved axis length
    n_ker, n_k_c = ker.shape

    # Algorithm works for odd kernel grid
    if n_ker % 2 != 1:
        raise ValueError("length of the convolution kernel should be odd.")

    # Assign half-length
    h_len = (n_ker - 1) // 2

    # Define each diagonal of the sparse convolution matrix
    diag_val, offset = [], []
    for i_ker, i_k_c in enumerate(range(-h_len, h_len+1)):

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
    """
    Return a convolution matrix
    Can return a sparse matrix (N_k_convolved, N_k)
    or a matrix in the compact form (N_ker, N_k_convolved).
    N_k is the length of the grid on which the convolution
    will be applied, N_k_convolved is the length of the
    grid after convolution and N_ker is the maximum length of
    the kernel. If the default sparse matrix option is chosen,
    the convolution can be apply on an array f | f = fct(grid)
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
        kernel = to_2d(kernel, grid, [a, b])
    else:
        # TODO add other options.
        pass

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
    """
    Define a gaussian convolution kernel at the nyquist
    sampling. For a given point on the grid x_i, the kernel
    is given by a gaussian with
    FWHM = n_sampling * (dx_(i-1) + dx_i) / 2.
    The FWHM is computed for each elements of the grid except
    the extremities (not defined). We can then generate FWHM as
    a function of thegrid and interpolate/extrapolate to get
    the kernel as a function of its position relative to the grid.
    """
    def __init__(self, grid, n_sampling=2, bounds_error=False,
                 fill_value="extrapolate", **kwargs):
        """
        Parameters
        ----------
        grid : 1d array
            Grid used to define the kernels
        n_sampling: int, optional
            sampling of the grid. Default is 2, so we assume that
            the grid is Nyquist sampled.
        bounds_error, fill_value and kwargs:
            `interp1d` kwargs used to get FWHM as a function of the grid.
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
        """
        Parameters
        ----------
        x: 1d array
            position where the kernel is evaluated
        x0: 1d array (same shape as x)
            position of the kernel center for each x.

        Returns
        -------
        Value of the gaussian kernel for each sets of (x, x0)
        """

        # Get the sigma of each gaussians
        sig = self.fct_sig(x0)

        return gaussians(x, x0, sig)


# ==============================================================================
# Code for doing Tikhonov regularisation.
# ==============================================================================


def finite_diff(x):
    """
    Returns the finite difference matrix operator based on x.
    Input:
        x: array-like
    Output:
        sparse matrix. When apply to x `diff_matrix.dot(x)`,
        the result is the same as np.diff(x)
    """
    n_x = len(x)

    # Build matrix
    diff_matrix = diags([-1.], shape=(n_x-1, n_x))
    diff_matrix += diags([1.], 1, shape=(n_x-1, n_x))

    return diff_matrix


def finite_second_d(grid):
    """
    Returns the second derivative operator based on grid
    Inputs:
    -------
    grid: 1d array-like
        grid where the second derivative will be compute.
    Ouputs:
    -------
    second_d: matrix
        Operator to compute the second derivative, so that
        f" = second_d.dot(f), where f is a function
        projected on `grid`.
    """

    # Finite difference operator
    d_matrix = finite_diff(grid)

    # Delta lambda
    d_grid = d_matrix.dot(grid)

    # First derivative operator
    first_d = diags(1./d_grid).dot(d_matrix)

    # Second derivative operator
    second_d = finite_diff(grid[:-1]).dot(first_d)

    # don't forget the delta labda
    second_d = diags(1./d_grid[:-1]).dot(second_d)

    return second_d


def finite_first_d(grid):
    """
    Returns the first derivative operator based on grid
    Inputs:
    -------
    grid: 1d array-like, optional
        grid where the first derivative will be compute.
    Ouputs:
    -------
    first_d: matrix
        Operator to compute the second derivative, so that
        f' = first_d.dot(f), where f is a function
        projected on `grid`.
    """

    # Finite difference operator
    d_matrix = finite_diff(grid)

    # Delta lambda
    d_grid = d_matrix.dot(grid)

    # First derivative operator
    first_d = diags(1./d_grid).dot(d_matrix)

    return first_d


def finite_zeroth_d(grid):
    """
    Gives the zeroth derivative operator on the function
    f(grid), so simply returns the identity matrix... XD
    """
    return identity(len(grid))


def get_nyquist_matrix(grid, integrate=True, n_sampling=2,
                       thresh=1e-5, **kwargs):
    """
    Get the tikhonov regularisation matrix based on
    a Nyquist convolution matrix (convolution with
    a kernel with a resolution given by the sampling
    of a grid). The Tikhonov matrix will be given by
    the difference of the nominal solution and
    the convolved solution.

    Parameters
    ----------
    grid: 1d-array
        Grid to project the kernel
    integrate: bool, optional
        If True, add integration weights to the tikhonov matrix, so
        when the squared norm is computed, the result is equivalent
        to the integral of the integrand squared.
    n_sampling: int, optional
        sampling of the grid. Default is 2, so we assume that
        the grid is Nyquist sampled.
    thresh: float, optional
        Used to define the maximum length of the kernel.
        Truncate when `kernel` < `thresh`
    kwargs:
        `interp1d` kwargs used to get FWHM as a function of the grid.
    """

    # Get nyquist kernel function
    ker = NyquistKer(grid, n_sampling=n_sampling, **kwargs)

    # Build convolution matrix
    conv_matrix = get_c_matrix(ker, grid, thresh=thresh)

    # Build tikhonov matrix
    t_mat = conv_matrix - identity(conv_matrix.shape[0])

    if integrate:
        # The grid may not be evenly spaced, so
        # add an integration weight
        d_grid = np.diff(grid)
        d_grid = np.concatenate([d_grid, [d_grid[-1]]])
        t_mat = diags(np.sqrt(d_grid)).dot(t_mat)

    return t_mat


def tikho_solve(a_mat, b_vec, t_mat=None, grid=None,
                verbose=True, factor=1.0, estimate=None, index=None):
    """
    Tikhonov solver to use as a function instead of a class.

    Parameters
    ----------
    a_mat: matrix-like object (2d)
        matrix A in the system to solve A.x = b
    b_vec: vector-like object (1d)
        vector b in the system to solve A.x = b
    t_mat: matrix-like object (2d), optional
        Tikhonov regularisation matrix to be applied on b_vec.
        Default is the default of the Tikhonov class. (Identity matrix)
    grid: array-like 1d, optional
        grid on which b-vec is projected. Used to compute derivative
    verbose: bool
        Print details or not
    factor: float, optional
        multiplicative constant of the regularisation matrix
    estimate: vector-like object (1d)
        Estimate oof the solution of the system.
    index: indexable, optional
        index of the valid row of the b_vec.

    Returns
    ------
    Solution of the system (1d array)
    """
    tikho = Tikhonov(a_mat, b_vec, t_mat=t_mat,
                     grid=grid, verbose=verbose, index=index)

    return tikho.solve(factor=factor, estimate=estimate)


class Tikhonov:
    """
    Tikhonov regularisation to solve the ill-condition problem:
    A.x = b, where A is accidently singular or close to singularity.
    Tikhonov regularisation adds a regularisation term in
    the equation and aim to minimize the equation:
    ||A.x - b||^2 + ||gamma.x||^2
    Where gamma is the Tikhonov regularisation matrix.
    """
    default_mat = {'zeroth': finite_zeroth_d,
                   'first': finite_first_d,
                   'second': finite_second_d}

    def __init__(self, a_mat, b_vec, t_mat=None,
                 grid=None, verbose=True, index=None):
        """
        Parameters
        ----------
        a_mat: matrix-like object (2d)
            matrix A in the system to solve A.x = b
        b_vec: vector-like object (1d)
            vector b in the system to solve A.x = b
        t_mat: matrix-like object (2d), optional
            Tikhonov regularisation matrix to be applied on b_vec.
            Default is the default of the Tikhonov class. (Identity matrix)
        grid: array-like 1d, optional
            grid on which b-vec is projected. Used to compute derivative
        verbose: bool
            Print details or not
        index: indexable, optional
            index of the valid row of the b_vec.
        """

        # b_vec will be passed to default_mat functions
        # if grid not given.
        if grid is None and t_mat is None:
            grid = b_vec

        # If string, search in the default Tikhonov matrix
        if t_mat is None:
            self.type = 'zeroth'
            t_mat = self.default_mat['zeroth'](grid)
        elif callable(t_mat):
            t_mat = t_mat(grid)
            self.type = 'custom'
        else:
            self.type = 'custom'

        # Take all indices by default
        if index is None:
            index = slice(None)

        self.a_mat = a_mat[index, :][:, index]
        self.b_vec = b_vec[index]
        self.t_mat = t_mat[index, :][:, index]
        self.index = index
        self.verbose = verbose
        self.test = None

        return

    def verbose_print(self, *args, **kwargs):
        """Print if verbose is True. Same as `print` function."""

        if self.verbose:
            print(*args, **kwargs)

        return

    def solve(self, factor=1.0, estimate=None):
        """
        Minimize the equation ||A.x - b||^2 + ||gamma.x||^2
        by solving (A_T.A + gamma_T.gamma).x = A_T.b
        gamma is the Tikhonov matrix multiplied by a scale factor

        Parameters
        ----------
        factor: float, optional
            multiplicative constant of the regularisation matrix
        estimate: vector-like object (1d)
            Estimate oof the solution of the system.

        Returns
        ------
        Solution of the system (1d array)
        """
        # Get needed attributes
        a_mat = self.a_mat
        b_vec = self.b_vec
        index = self.index

        # Matrix gamma (with scale factor)
        gamma = factor * self.t_mat

        # Build system
        gamma_2 = (gamma.T).dot(gamma)  # Gamma square
        matrix = a_mat.T.dot(a_mat) + gamma_2
        result = (a_mat.T).dot(b_vec.T)

        # Include solution estimate if given
        if estimate is not None:
            result += gamma_2.dot(estimate[index].T)

        # Solve
        return spsolve(matrix, result)

    def test_factors(self, factors, estimate=None):
        """
        test multiple factors

        Parameters
        ----------
        factors: 1d array-like
            factors to test
        estimate: array like
            estimate of the solution

        Returns
        ------
        dictionnary of test results
        """

        self.verbose_print('Testing factors...')

        # Get relevant attributes
        b_vec = self.b_vec
        a_mat = self.a_mat
        t_mat = self.t_mat

        # Init outputs
        sln, err, reg = [], [], []

        # Test all factors
        for i_fac, factor in enumerate(factors):

            # Save solution
            sln.append(self.solve(factor, estimate))

            # Save error A.x - b
            err.append(a_mat.dot(sln[-1]) - b_vec)

            # Save regulatisation term
            reg.append(t_mat.dot(sln[-1]))

            # Print
            message = '{}/{}'.format(i_fac, len(factors))
            self.verbose_print(message, end='\r')

        # Final print
        message = '{}/{}'.format(i_fac + 1, len(factors))
        self.verbose_print(message)

        # Convert to arrays
        sln = np.array(sln)
        err = np.array(err)
        reg = np.array(reg)

        # Save in a dictionnary
        self.test = {'factors': factors,
                     'solution': sln,
                     'error': err,
                     'reg': reg}

        return self.test

    def _check_plot_inputs(self, fig, ax, label, factors, test):
        """
        Method to manage inputs for plots methods.
        """

        # Use ax or fig if given. Else, init the figure
        if (fig is None) and (ax is None):
            fig, ax = plt.subplots(1, 1, sharex=True)
        elif ax is None:
            ax = fig.subplots(1, 1, sharex=True)

        # Use the type of regularisation as label if None is given
        if label is None:
            label = self.type

        if test is None:

            # Run tests with `factors` if not done already.
            if self.test is None:
                self.test_factors(factors)

            test = self.test

        return fig, ax, label, test

    def error_plot(self, fig=None, ax=None, factors=None,
                   label=None, test=None, test_key=None, y_val=None):
        """
        Plot error as a function of factors

        Parameters
        ----------
        fig: matplotlib figure, optional
            Figure to use for plot
            If not given and ax is None, new figure is initiated
        ax: matplotlib axis, optional
            axis to use for plot. If not given, a new axis is initiated.
        factors: 1d array-like
            factors to test
        label: str, optional
            label too put in legend
        test: dictionnary, optional
            dictionnary of tests (output of Tikhonov.test_factors)
        test_key: str, optional
            which test result to plot. If not specified,
            the euclidian norm of the 'error' key will be used.
        y_val: array-like, optional
            y values to plot. Same length as factors.

        Returns
        ------
        fig, ax
        """

        # Manage method's inputs
        args = (fig, ax, label, factors, test)
        fig, ax, label, test = self._check_plot_inputs(*args)

        # What y value do we plot?
        if y_val is None:

            # Use tests to plot y_val
            if test_key is None:

                # Default is euclidian norm of error.
                # Similar to the chi^2.
                y_val = (test['error']**2).sum(axis=-1)
            else:
                y_val = test[test_key]

        # Plot
        ax.loglog(test['factors'], y_val, label=label)

        # Mark minimum value
        i_min = np.argmin(y_val)
        min_coord = test['factors'][i_min], y_val[i_min]
        ax.scatter(*min_coord, marker="x")
        text = '{:2.1e}'.format(min_coord[0])
        ax.text(*min_coord, text, va="top", ha="center")

        # Show legend
        ax.legend()

        # Labels
        ax.set_xlabel("Scale factor")
        ylabel = r'System error '
        ylabel += r'$\left(||\mathbf{Ax-b}||^2_2\right)$'
        ax.set_ylabel(ylabel)

        return fig, ax

    def l_plot(self, fig=None, ax=None, factors=None, label=None,
               test=None, text_label=True, factor_norm=False):
        """
        make an 'l plot'

        Parameters
        ----------
        fig: matplotlib figure, optional
            Figure to use for plot
            If not given and ax is None, new figure is initiated
        ax: matplotlib axis, optional
            axis to use for plot. If not given, a new axis is initiated.
        factors: 1d array-like
            factors to test
        label: str, optional
            label too put in legend
        test: dictionnary, optional
            dictionnary of tests (output of Tikhonov.test_factors)
        text_label: bool, optional
            Add label of the factor value to each points in the plot.

        Returns
        ------
        fig, ax
        """

        # Manage method's inputs
        args = (fig, ax, label, factors, test)
        fig, ax, label, test = self._check_plot_inputs(*args)

        # Compute euclidian norm of error (||A.x - b||).
        # Similar to the chi^2.
        err_norm = (test['error']**2).sum(axis=-1)

        # Compute norm of regularisation term
        reg_norm = (test['reg']**2).sum(axis=-1)

        # Factors
        if factor_norm:
            reg_norm *= test['factors']**2

        # Plot
        ax.loglog(err_norm, reg_norm, '.:', label=label)

        # Add factor values as text
        if text_label:
            for f, x, y in zip(test['factors'], err_norm, reg_norm):
                plt.text(x, y, "{:2.1e}".format(f), va="center", ha="right")

        # Legend
        ax.legend()

        # Labels
        xlabel = r'$\left(||\mathbf{Ax-b}||^2_2\right)$'
        ylabel = r'$\left(||\mathbf{\Gamma.x}||^2_2\right)$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return fig, ax


def main():

    return


if __name__ == '__main__':
    main()
