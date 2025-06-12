import numpy as np
import pytest
import logging
from scipy.integrate import trapezoid

from jwst.tests.helpers import LogWatcher
from jwst.extract_1d.soss_extract import atoca_utils as au


def test_arange_2d():
    starts = np.array([3, 4, 5])
    stops = np.ones(starts.shape) * 7
    out = au.arange_2d(starts, stops)

    bad = -1
    expected_out = np.array([[3, 4, 5, 6], [4, 5, 6, bad], [5, 6, bad, bad]])
    assert np.allclose(out, expected_out)

    # test bad input catches
    starts_wrong_shape = starts[1:]
    with pytest.raises(ValueError):
        au.arange_2d(starts_wrong_shape, stops)

    stops_too_small = np.copy(stops)
    stops_too_small[2] = 4
    with pytest.raises(ValueError):
        au.arange_2d(starts, stops_too_small)


@pytest.mark.parametrize("dispersion_axis", [0, 1])
def test_get_wv_map_bounds(wave_map, dispersion_axis):
    """
    top is the low-wavelength end, bottom is high-wavelength end
    """
    wave_map = wave_map[0].copy()
    wave_map[1, 3] = -1  # test skip of bad value
    wavelengths = wave_map[0]

    if dispersion_axis == 0:
        wave_flip = wave_map.T
    else:
        wave_flip = wave_map
    wave_top, wave_bottom = au._get_wv_map_bounds(wave_flip, dispersion_axis=dispersion_axis)

    # flip the results back so we can reuse the same tests
    if dispersion_axis == 0:
        wave_top = wave_top.T
        wave_bottom = wave_bottom.T

    diff = (wavelengths[1:] - wavelengths[:-1]) / 2
    diff_lower = np.insert(diff, 0, diff[0])
    diff_upper = np.append(diff, diff[-1])
    wave_top_expected = wavelengths - diff_lower
    wave_bottom_expected = wavelengths + diff_upper

    # basic test
    assert wave_top.shape == wave_bottom.shape == (wave_map.shape[0],) + wavelengths.shape
    assert np.allclose(wave_top[0], wave_top_expected)
    assert np.allclose(wave_bottom[0], wave_bottom_expected)

    # test skip bad pixel
    assert wave_top[1, 3] == 0
    assert wave_bottom[1, 3] == 0

    # test bad input error raises
    with pytest.raises(ValueError):
        au._get_wv_map_bounds(wave_flip, dispersion_axis=2)


@pytest.mark.parametrize("dispersion_axis", [0, 1])
def test_get_wave_p_or_m(wave_map, dispersion_axis):
    """
    Check that the plus and minus side is correctly identified
    for strictly ascending and strictly descending wavelengths.
    """
    wave_map = wave_map[0].copy()
    wave_reverse = np.fliplr(wave_map)
    if dispersion_axis == 0:
        wave_flip = wave_map.T
        wave_reverse = wave_reverse.T
    else:
        wave_flip = wave_map

    wave_p_0, wave_m_0 = au.get_wave_p_or_m(wave_flip, dispersion_axis=dispersion_axis)
    wave_p_1, wave_m_1 = au.get_wave_p_or_m(wave_reverse, dispersion_axis=dispersion_axis)

    if dispersion_axis == 0:
        wave_p_0 = wave_p_0.T
        wave_m_0 = wave_m_0.T
        wave_p_1 = wave_p_1.T
        wave_m_1 = wave_m_1.T
    assert np.all(wave_p_0 >= wave_m_0)
    assert np.allclose(wave_p_0, np.fliplr(wave_p_1))
    assert np.allclose(wave_m_0, np.fliplr(wave_m_1))


def test_get_wave_p_or_m_not_ascending(wave_map):
    wave_map = wave_map[0].copy()
    with pytest.raises(ValueError):
        wave_map[0, 5] = 2  # make it not strictly ascending
        au.get_wave_p_or_m(wave_map, dispersion_axis=1)


FIBONACCI = np.array([1, 1, 2, 3, 5, 8, 13, 21, 35], dtype=float)


@pytest.mark.parametrize("n_os", [1, 5])
def test_oversample_grid(n_os):
    oversample = au.oversample_grid(FIBONACCI, n_os)

    # oversample_grid is supposed to remove any duplicates, and there is a duplicate
    # in FIBONACCI. So the output should be 4 times the size of FIBONACCI minus 1
    assert oversample.size == n_os * (FIBONACCI.size - 1) - (n_os - 1)
    assert oversample.min() == FIBONACCI.min()
    assert oversample.max() == FIBONACCI.max()

    # test whether np.interp could have been used instead
    grid = np.arange(0, FIBONACCI.size, 1 / n_os)
    wls = np.unique(np.interp(grid, np.arange(FIBONACCI.size), FIBONACCI))
    assert np.allclose(oversample, wls)


@pytest.mark.parametrize("os_factor", [1, 2, 5])
def test_oversample_irregular(os_factor):
    """Test oversampling to a grid with irregular spacing"""
    # oversampling function removes duplicates,
    # this is tested in previous test, and just complicates counting for this test
    # for FIBONACCI, unique is just removing zeroth element
    fib_unq = np.unique(FIBONACCI)
    n_os = np.ones((fib_unq.size - 1,), dtype=int)
    n_os[2:5] = os_factor
    n_os[3] = os_factor * 2
    # this gives n_os = [1 1 2 4 2] for os_factor = 2

    oversample = au.oversample_grid(fib_unq, n_os)

    # test no oversampling was done on the elements where not requested
    assert np.allclose(oversample[0:2], fib_unq[0:2])
    assert np.allclose(oversample[-1:], fib_unq[-1:])

    # test output shape.
    assert oversample.size == np.sum(n_os) + 1

    # test that this could have been done easily with np.interp
    intervals = 1 / n_os
    intervals = np.insert(np.repeat(intervals, n_os), 0, 0)
    grid = np.cumsum(intervals)
    wls = np.interp(grid, np.arange(fib_unq.size), fib_unq)
    assert wls.size == oversample.size
    assert np.allclose(oversample, wls)

    # test that n_os shape must match input shape - 1
    with pytest.raises(ValueError):
        au.oversample_grid(fib_unq, n_os[:-1])


WAVELENGTHS = np.linspace(1.5, 3.0, 50) + np.sin(np.linspace(0, np.pi / 2, 50))


@pytest.mark.parametrize("wave_range", [(2.1, 3.9), (1.8, 4.5)])
def test_extrapolate_grid(wave_range):
    extrapolated = au._extrapolate_grid(WAVELENGTHS, wave_range, 1)

    assert extrapolated.max() > wave_range[1]
    assert extrapolated.min() < wave_range[0]
    assert np.all(extrapolated[1:] >= extrapolated[:-1])

    # if interpolation not needed on either side, should return the original
    if wave_range == (2.1, 3.9):
        assert extrapolated is WAVELENGTHS


def test_extrapolate_catch_failed_converge():
    # give wavelengths some non-linearity
    wave_range = WAVELENGTHS.min(), WAVELENGTHS.max() + 4.0
    with pytest.raises(RuntimeError):
        au._extrapolate_grid(WAVELENGTHS, wave_range, 1)


def test_extrapolate_bad_inputs():
    with pytest.raises(ValueError):
        au._extrapolate_grid(WAVELENGTHS, (2.9, 2.1))
    with pytest.raises(ValueError):
        au._extrapolate_grid(WAVELENGTHS, (4.1, 4.2))


def test_grid_from_map(wave_map, trace_profile):
    """
    Covers expected behavior of grid_from_map_with_extrapolation,
    including coverage of a previous bug
    where bad wavelengths were not being ignored properly.
    """

    wave_map = wave_map[0].copy()
    wavelengths = wave_map[0][::-1]
    trace_profile = trace_profile[0].copy()
    wave_grid = au.grid_from_map_with_extrapolation(wave_map, trace_profile, wave_range=None)

    assert np.allclose(wave_grid, wavelengths)

    # test custom wave_range
    wave_range = [wavelengths[2], wavelengths[-2] + 0.01]
    wave_grid = au.grid_from_map_with_extrapolation(wave_map, trace_profile, wave_range=wave_range)
    assert np.allclose(wave_grid, wavelengths[2:-1])

    # test custom wave_range with extrapolation
    wave_range = [wavelengths[2], wavelengths[-1] + 1]
    wave_grid = au.grid_from_map_with_extrapolation(wave_map, trace_profile, wave_range=wave_range)
    assert len(wave_grid) > len(wavelengths[2:])
    n_inside = wavelengths[2:].size
    assert np.allclose(wave_grid[:n_inside], wavelengths[2:])

    with pytest.raises(ValueError):
        au.grid_from_map_with_extrapolation(wave_map, trace_profile, wave_range=[0.1, 0.2])


def xsinx(x):
    return x * np.sin(x)


def test_estim_integration_error():
    """
    Use as truth the x sin(x) from 0 to pi, has an analytic solution == pi.
    """

    n = 11
    grid = np.linspace(0, np.pi, n)
    err, rel_err = au._estim_integration_err(grid, xsinx)

    assert len(rel_err) == n - 1
    assert np.all(rel_err >= 0)
    assert np.all(rel_err < 1)


@pytest.mark.parametrize("max_iter, rtol", [(1, 1e-3), (10, 1e-9), (10, 1e-3), (1, 1e-9)])
def test_adapt_grid(max_iter, rtol):
    """
    Use as truth the x sin(x) from 0 to pi, has an analytic solution == pi.
    """

    input_grid = np.linspace(0, np.pi, 11)
    input_grid_diff = input_grid[1] - input_grid[0]
    max_grid_size = 100
    grid, is_converged = au._adapt_grid(
        input_grid, xsinx, max_grid_size, max_iter=max_iter, rtol=rtol
    )

    # ensure grid respects max_grid_size and max_iter in all cases
    assert len(grid) <= max_grid_size
    grid_diff = grid[1:] - grid[:-1]
    assert np.min(grid_diff) >= input_grid_diff / (2**max_iter)

    numerical_integral = trapezoid(xsinx(grid), grid)

    # ensure this converges for at least one of our test cases
    if max_iter == 10 and rtol == 1e-3:
        assert is_converged

    if is_converged:
        # test error of the answer is smaller than rtol
        assert np.isclose(numerical_integral, np.pi, rtol=rtol)
        # test that success was a stop condition
        assert len(grid) < max_grid_size

    # test stop conditions
    elif max_iter == 10:
        # ensure hitting max_grid_size returns an array of exactly length max_grid_size
        assert len(grid) == max_grid_size
    elif max_iter == 1:
        # ensure hitting max_iter can stop iteration before max_grid_size reached
        assert len(grid) <= 2 * len(input_grid)


def test_adapt_grid_bad_inputs():
    with pytest.raises(ValueError):
        # input grid larger than max_grid_size
        au._adapt_grid(np.array([1, 2, 3]), xsinx, 2)


def test_trim_grids():
    grid_range = (-3, 3)
    grid0 = np.linspace(-3, 0, 4)  # kept entirely.
    grid1 = np.linspace(
        -3, 0, 16
    )  # removed entirely. Finer spacing doesn't matter, preceded by grid0
    grid2 = np.linspace(0, 3, 5)  # kept from 0 to 3
    grid3 = np.linspace(
        -4, 4, 5
    )  # removed entirely. Outside of grid_range and the rest is superseded

    all_grids = [grid0, grid1, grid2, grid3]
    trimmed_grids = au._trim_grids(all_grids, grid_range)

    assert len(trimmed_grids) == len(all_grids)
    assert trimmed_grids[0].size == grid0.size
    assert trimmed_grids[1].size == 0
    assert trimmed_grids[2].size == grid2.size
    assert trimmed_grids[3].size == 0


def test_make_combined_adaptive_grid():
    """see also tests of _adapt_grid and _trim_grids for more detailed tests"""

    grid_range = (0, np.pi)
    grid0 = np.linspace(0, np.pi / 2, 6)  # kept entirely.
    grid1 = np.linspace(
        0, np.pi / 2, 15
    )  # removed entirely. Finer spacing doesn't matter, preceded by grid0
    grid2 = np.linspace(np.pi / 2, np.pi, 11)  # kept from pi/2 to pi

    # purposely make same lower index for grid2 as upper index for grid0 to test uniqueness of output

    all_grids = [grid0, grid1, grid2]
    all_estimate = [xsinx, xsinx, xsinx]

    rtol = 1e-3
    combined_grid = au.make_combined_adaptive_grid(
        all_grids, all_estimate, grid_range, max_iter=10, rtol=rtol, max_total_size=100
    )

    numerical_integral = trapezoid(xsinx(combined_grid), combined_grid)

    assert np.unique(combined_grid).size == combined_grid.size
    assert np.isclose(numerical_integral, np.pi, rtol=rtol)


def test_make_combined_adaptive_grid_maxsize(monkeypatch):
    """
    Test that max_total_size is respected.

    It's a bit more complicated than just max_total_size is the size of the output,
    because the code allows the combined grid to add on grid2 at its native resolution
    if grid0 fails to converge and reaches the max_total_size.
    """
    # make a log watcher
    watcher = LogWatcher("Precision cannot be guaranteed: max grid size")
    monkeypatch.setattr(
        logging.getLogger("jwst.extract_1d.soss_extract.atoca_utils"), "warning", watcher
    )

    grid_range = (0, np.pi)
    grid0 = np.linspace(0, np.pi / 2, 6)  # kept entirely.
    grid2 = np.linspace(np.pi / 2 + 0.001, np.pi, 11)  # kept from pi/2 to pi
    max_grid_size = 100
    all_grids = [grid0, grid2]
    all_estimate = [xsinx, xsinx]
    combined_grid = au.make_combined_adaptive_grid(
        all_grids, all_estimate, grid_range, max_iter=10, rtol=1e-4, max_total_size=max_grid_size
    )
    assert combined_grid.size == max_grid_size + grid2.size
    watcher.assert_seen()


def test_throughput_soss():
    wavelengths = np.linspace(2, 5, 10)
    throughputs = np.ones_like(wavelengths)
    interpolator = au.throughput_soss(wavelengths, throughputs)

    # test that it returns 1 for all wavelengths inside range
    interp = interpolator(wavelengths)
    assert np.allclose(interp[1:-1], throughputs[1:-1])
    assert interp[0] == 0
    assert interp[-1] == 0

    # test that it returns 0 for all wavelengths outside range
    wavelengths_outside = np.linspace(1, 1.5, 5)
    interp = interpolator(wavelengths_outside)
    assert np.all(interp == 0)

    # test ValueError raise for shape mismatch
    with pytest.raises(ValueError):
        au.throughput_soss(wavelengths, throughputs[:-1])


def test_webb_kernel(webb_kernels, wave_map):
    wave_trace = wave_map[0][0]
    min_trace, max_trace = np.min(wave_trace), np.max(wave_trace)
    kern = webb_kernels[0]

    # basic ensure that the input is stored and shapes
    assert kern.wave_kernels.shape == kern.kernels.shape

    # test that pixels and wave_kernels are both monotonic
    assert np.all(np.diff(kern.pixels) > 0)
    assert np.all(np.diff(kern.wave_kernels) > 0)

    # test that pixels is mirrored around the center and has zero at center
    assert np.allclose(kern.pixels + kern.pixels[::-1], 0)
    assert kern.pixels[kern.pixels.size // 2] == 0

    # test that wave_center has same shape as wavelength axis of wave_kernel
    # but contains values that are in wave_trace
    assert kern.wave_center.size == kern.wave_kernels.shape[1]
    assert all(np.isin(kern.wave_center, wave_trace))

    # test min value
    assert kern.min_value > 0
    assert np.isin(kern.min_value, kern.kernels)
    assert isinstance(kern.min_value, float)

    # test the polynomial fit has the proper shape. hard-coded to a first-order, i.e., linear fit
    # since the throughput is constant in wavelength, the slopes should be close to zero
    # and the y-intercepts should be close to kern.wave_center
    # especially with so few points. just go with 10 percent, should catch egregious changes
    assert kern.poly.shape == (kern.wave_kernels.shape[1], 2)
    assert np.allclose(kern.poly[:, 0], 0, atol=1e-1)
    assert np.allclose(kern.poly[:, 1], kern.wave_center, atol=1e-1)

    # test interpolation function, which takes in a pixel and a wavelength and returns a throughput
    # this should return the triangle function at all wavelengths and zero outside range
    pix_half = kern.n_pix // 2
    wl_test = np.linspace(min_trace, max_trace, 10)
    pixels_test = np.array([-pix_half - 1, 0, pix_half, pix_half + 1])

    data_in = kern.kernels[:, 0]
    m = kern.min_value
    expected = np.array([m, np.max(data_in), m, m])

    interp = kern.f_ker(pixels_test, wl_test)
    assert interp.shape == (pixels_test.size, wl_test.size)
    diff = interp[:, 1:] - interp[:, :-1]
    assert np.allclose(diff, 0)
    assert np.allclose(interp[:, 0], expected, rtol=1e-3)

    # call the kernel object directly
    # this takes a wavelength and a central wavelength of the kernel,
    # then converts to pixels to use self.f_ker internally
    kern_val = kern(wl_test, wl_test)
    assert kern(wl_test, wl_test).ndim == 1
    # ignore edge effects
    assert np.allclose(kern_val[1:-1], np.max(data_in))

    # both inputs need to be same shape
    with pytest.raises(ValueError):
        kern(wl_test, wl_test[:-1])


def test_finite_first_diff():
    wave_grid = np.linspace(0, 2 * np.pi, 100)
    test_0 = np.ones_like(wave_grid)
    test_sin = np.sin(wave_grid)

    first_d = au.finite_first_d(wave_grid)
    assert first_d.size == (wave_grid.size - 1) * 2

    # test trivial example returning zeros for constant
    f0 = first_d.dot(test_0)
    assert np.allclose(f0, 0)

    # test derivative of sin returns cos
    wave_between = (wave_grid[1:] + wave_grid[:-1]) / 2
    f_sin = first_d.dot(test_sin)
    assert np.allclose(f_sin, np.cos(wave_between), atol=1e-3)


def test_get_c_matrix(kernels_unity, webb_kernels, wave_grid):
    """See also test_fct_to_array and test_sparse_c for more detailed tests
    of functions called by this one"""

    # only need to test one order
    kern = webb_kernels[0]
    matrix = au.get_c_matrix(kern, wave_grid, i_bounds=None)

    # ensure proper shape
    assert matrix.shape == (wave_grid.size, wave_grid.size)
    assert matrix.dtype == np.float64

    # ensure normalized
    assert np.isclose(matrix.sum(), matrix.shape[0])

    # test where input kernel is a 2-D array instead of callable
    i_bounds = [0, len(wave_grid)]
    kern_array = au._fct_to_array(kern, wave_grid, i_bounds, 1e-5)
    matrix_from_array = au.get_c_matrix(kern_array, wave_grid, i_bounds=i_bounds)
    assert np.allclose(matrix.toarray(), matrix_from_array.toarray())

    # test where input kernel is size 1
    kern_unity = kernels_unity[0]
    matrix_from_unity = au.get_c_matrix(kern_unity, wave_grid, i_bounds=i_bounds)
    assert matrix_from_unity.shape == (wave_grid.size, wave_grid.size)

    # test where i_bounds is not None
    i_bounds = [10, wave_grid.size - 10]
    matrix_ibnds = au.get_c_matrix(kern, wave_grid, i_bounds=i_bounds)
    expected_shape = (wave_grid[i_bounds[0] : i_bounds[1]].size, wave_grid.size)
    assert matrix_ibnds.shape == expected_shape

    # Test invalid kernel input (wrong dimensions)
    with pytest.raises(ValueError):
        kern_array_bad = kern_array[np.newaxis, ...]
        au.get_c_matrix(kern_array_bad, wave_grid, i_bounds=i_bounds)

    # Test invalid kernel input (odd shape)
    with pytest.raises(ValueError):
        kern_array_bad = kern_array[1:, 1:]
        au.get_c_matrix(kern_array_bad, wave_grid, i_bounds=i_bounds)
