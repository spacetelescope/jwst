import pytest
import numpy as np
from functools import partial
from jwst.extract_1d.soss_extract import atoca
from jwst.extract_1d.soss_extract import atoca_utils as au


"""
Tests for atoca.py and atoca_utils.py

Create a miniature, slightly simplified model of the SOSS detector/optics.
Use those to instantiate an extraction engine, use the engine to create mock data
from a known input spectrum, and then check if the engine can retrieve that spectrum
from the data.

The model has the following features:
- Factor-of-10 smaller along each dimension
- Similar wavelengths for each of the two orders
- Partially overlapping traces for the two orders
- Randomly-selected bad pixels in the data
- Wave grid of size ~100 with varying resolution
- Triangle function throughput for each spectral order
- Kernel set to unity (for now)
"""


DATA_SHAPE = (25,200)
WAVE_BNDS_O1 = [2.8, 0.8]
WAVE_BNDS_O2 = [1.4, 0.5]
WAVE_BNDS_GRID = [0.7, 2.7]
ORDER1_SCALING = 20.0
ORDER2_SCALING = 2.0
SPECTRAL_SLOPE = 2


@pytest.fixture(scope="module")
def wave_map():
    """Trying for a roughly factor-of-10 shrink on each axis
    but otherwise relatively faithful reproduction of real detector behavior"""
    wave_ord1 = np.linspace(WAVE_BNDS_O1[0], WAVE_BNDS_O1[1], DATA_SHAPE[1])
    wave_ord1 = np.ones(DATA_SHAPE)*wave_ord1[np.newaxis, :]

    wave_ord2 = np.linspace(WAVE_BNDS_O2[0], WAVE_BNDS_O2[1], DATA_SHAPE[1])
    wave_ord2 = np.ones(DATA_SHAPE)*wave_ord2[np.newaxis,:]
    # add a small region of zeros to mimic what is input to the step from ref files
    wave_ord2[:,190:] = 0.0

    return [wave_ord1, wave_ord2]


@pytest.fixture(scope="module")
def trace_profile(wave_map):
    """
    order 2 is partially on top of, partially not on top of order 1
    give order 2 some slope to simulate that
    """
    # order 1
    DATA_SHAPE = wave_map[0].shape
    ord1 = np.zeros((DATA_SHAPE[0]))
    ord1[3:9] = 0.2
    ord1[2] = 0.1
    ord1[9] = 0.1
    profile_ord1 = np.ones(DATA_SHAPE)*ord1[:, np.newaxis]

    # order 2
    yy, xx = np.meshgrid(np.arange(DATA_SHAPE[0]), np.arange(DATA_SHAPE[1]))
    yy = yy.astype(np.float32) - xx.astype(np.float32)*0.08
    yy = yy.T
    
    profile_ord2 = np.zeros_like(yy)
    full = (yy >= 3) & (yy < 9)
    half0 = (yy >= 9) & (yy < 11)
    half1 = (yy >= 1) & (yy < 3)
    profile_ord2[full] = 0.2
    profile_ord2[half0] = 0.1
    profile_ord2[half1] = 0.1

    return [profile_ord1, profile_ord2]


@pytest.fixture(scope="module")
def wave_grid():
    """wave_grid has smaller spacings in some places than others
    and is not backwards order like the wave map
    Two duplicates are in there on purpose for testing"""
    lo0 = np.linspace(WAVE_BNDS_GRID[0], 1.2, 16)
    hi = np.linspace(1.2, 1.7, 46)
    lo2 = np.linspace(1.7, WAVE_BNDS_GRID[1], 31)
    return np.concatenate([lo0, hi, lo2])


@pytest.fixture(scope="module")
def throughput():
    """make a triangle function for each order but with different peak wavelength
    """

    def filter_function(wl, wl_max):
        """Set free parameters to roughly mimic throughput functions on main"""
        maxthru = 0.4
        thresh = 0.01
        scaling = 0.3
        dist = np.abs(wl - wl_max)
        thru = maxthru - dist*scaling
        thru[thru<thresh] = thresh
        return thru
    
    thru_o1 = partial(filter_function, wl_max=1.7)
    thru_o2 = partial(filter_function, wl_max=0.7)

    return [thru_o1, thru_o2]

@pytest.fixture(scope="module")
def kernels_unity():
    """For now, just return unity"""
    return [np.array([1.,]), np.array([1.,])]


@pytest.fixture(scope="module")
def webb_kernels(wave_map):
    """
    Toy model of the JWST kernels.
    Let the kernel be a triangle function along the pixel position axis,
    peaking at the center,
    and independent of wavelength.

    Same for both orders except the wavelengths and wave_trace are different.
    """
    n_os = 5
    n_pix = 15 # full pixel width of kernel
    n_wave = 10
    peakiness = 4
    min_val = 1
    kernel_width = n_os*n_pix - (n_os - 1)
    ctr_idx = kernel_width//2
    
    kernels = []
    for order in [0,1]:
        # set up wavelength grid over which kernel is defined
        wave_trace = wave_map[order][0]
        wave_range = (np.min(wave_trace), np.max(wave_trace))
        wavelengths = np.linspace(*wave_range, n_wave)
        wave_kernel = np.ones((kernel_width, wavelengths.size), dtype=float)*wavelengths[None,:]

        # model kernel as simply a triangle function that peaks at the center
        triangle_function = ctr_idx - peakiness*np.abs(ctr_idx - np.arange(0, kernel_width))
        triangle_function[triangle_function<=min_val] = min_val
        kernel = np.ones((kernel_width, wavelengths.size), dtype=float)*triangle_function[:,None]
        kernel/=np.sum(kernel)
        
        kernels.append(au.WebbKernel(wave_kernel, kernel, wave_trace, n_pix))

    return kernels


@pytest.fixture(scope="module")
def mask_trace_profile(trace_profile):
    """Masks set to unity where bad"""

    def mask_from_trace(trace_in, cut_low=None, cut_hi=None):
        trace = trace_in.copy()
        trace[trace_in<=0] = 1
        trace[trace_in>0] = 0
        if cut_low is not None:
            trace[:,:cut_low] = 1
        if cut_hi is not None:
            trace[:,cut_hi:] = 1
        return trace.astype(bool)
    
    trace_o1 = mask_from_trace(trace_profile[0], cut_low=0, cut_hi=199)
    trace_o2 = mask_from_trace(trace_profile[1], cut_low=0, cut_hi=175)
    return [trace_o1, trace_o2]


@pytest.fixture(scope="module")
def detector_mask(wave_map):
    """Add a few random bad pixels"""
    rng = np.random.default_rng(42)
    mask = np.zeros(DATA_SHAPE, dtype=bool)
    bad = rng.choice(mask.size, 100)
    bad = np.unravel_index(bad, DATA_SHAPE)
    mask[bad] = 1
    return mask


@pytest.fixture(scope="module")
def engine(wave_map,
           trace_profile,
           throughput,
           kernels_unity,
           wave_grid,
           mask_trace_profile,
           detector_mask,
):
    return atoca.ExtractionEngine(wave_map,
                                    trace_profile,
                                    throughput,
                                    kernels_unity,
                                    wave_grid,
                                    mask_trace_profile,
                                    global_mask=detector_mask)


def test_arange_2d():

    starts = np.array([3,4,5])
    stops = np.ones(starts.shape)*7
    out = au.arange_2d(starts, stops)

    bad = -1
    expected_out = np.array([
        [3,4,5,6],
        [4,5,6,bad],
        [5,6,bad,bad]
    ])
    assert np.allclose(out, expected_out)

    # test bad input catches
    starts_wrong_shape = starts[1:]
    with pytest.raises(ValueError):
        au.arange_2d(starts_wrong_shape, stops)
    
    stops_too_small = np.copy(stops)
    stops_too_small[2] = 4
    with pytest.raises(ValueError):
        au.arange_2d(starts, stops_too_small)


@pytest.mark.parametrize("dispersion_axis", [0,1])
def test_get_wv_map_bounds(wave_map, dispersion_axis):
    """
    top is the low-wavelength end, bottom is high-wavelength end
    """
    wave_map = wave_map[0].copy()
    wave_map[1,3] = -1 #test skip of bad value
    wavelengths = wave_map[0]

    if dispersion_axis == 0:
        wave_flip = wave_map.T
    else:
        wave_flip = wave_map
    wave_top, wave_bottom = au._get_wv_map_bounds(wave_flip, dispersion_axis=dispersion_axis)

    # flip the results back so we can re-use the same tests
    if dispersion_axis == 0:
        wave_top = wave_top.T
        wave_bottom = wave_bottom.T

    diff = (wavelengths[1:]-wavelengths[:-1])/2
    diff_lower = np.insert(diff,0,diff[0])
    diff_upper = np.append(diff,diff[-1])
    wave_top_expected = wavelengths-diff_lower
    wave_bottom_expected = wavelengths+diff_upper

    # basic test
    assert wave_top.shape == wave_bottom.shape == (wave_map.shape[0],)+wavelengths.shape
    assert np.allclose(wave_top[0], wave_top_expected)
    assert np.allclose(wave_bottom[0], wave_bottom_expected)

    # test skip bad pixel
    assert wave_top[1,3] == 0
    assert wave_bottom[1,3] == 0

    # test bad input error raises
    with pytest.raises(ValueError):
        au._get_wv_map_bounds(wave_flip, dispersion_axis=2)


@pytest.mark.parametrize("dispersion_axis", [0,1])
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

    if dispersion_axis==0:
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
        wave_map[0,5] = 2 # make it not strictly ascending
        au.get_wave_p_or_m(wave_map, dispersion_axis=1)


FIBONACCI = np.array([1,1,2,3,5,8,13,21,35], dtype=float)
@pytest.mark.parametrize("n_os", [1,5])
def test_oversample_grid(n_os):

    oversample = au.oversample_grid(FIBONACCI, n_os)

    # oversample_grid is supposed to remove any duplicates, and there is a duplicate
    # in FIBONACCI. So the output should be 4 times the size of FIBONACCI minus 1
    assert oversample.size == n_os*(FIBONACCI.size - 1) - (n_os-1)
    assert oversample.min() == FIBONACCI.min()
    assert oversample.max() == FIBONACCI.max()

    # test whether np.interp could have been used instead
    grid = np.arange(0, FIBONACCI.size, 1/n_os)
    wls = np.unique(np.interp(grid, np.arange(FIBONACCI.size), FIBONACCI))
    assert np.allclose(oversample, wls)


@pytest.mark.parametrize("os_factor", [1,2,5])
def test_oversample_irregular(os_factor):
    """Test oversampling to a grid with irregular spacing"""
    # oversampling function removes duplicates,
    # this is tested in previous test, and just complicates counting for this test
    # for FIBONACCI, unique is just removing zeroth element
    fib_unq = np.unique(FIBONACCI) 
    n_os = np.ones((fib_unq.size-1,), dtype=int)
    n_os[2:5] = os_factor
    n_os[3] = os_factor*2
    # this gives n_os = [1 1 2 4 2] for os_factor = 2

    oversample = au.oversample_grid(fib_unq, n_os)

    # test no oversampling was done on the elements where not requested
    assert np.allclose(oversample[0:2], fib_unq[0:2])
    assert np.allclose(oversample[-1:], fib_unq[-1:])

    # test output shape.
    assert oversample.size == np.sum(n_os)+1

    # test that this could have been done easily with np.interp
    intervals = 1/n_os
    intervals = np.insert(np.repeat(intervals, n_os),0,0)
    grid = np.cumsum(intervals)
    wls = np.interp(grid, np.arange(fib_unq.size), fib_unq)
    assert wls.size == oversample.size
    assert np.allclose(oversample, wls)

    # test that n_os shape must match input shape - 1
    with pytest.raises(ValueError):
        au.oversample_grid(fib_unq, n_os[:-1])


WAVELENGTHS = np.linspace(1.5, 3.0, 50) + np.sin(np.linspace(0, np.pi/2, 50))
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
    wave_range = WAVELENGTHS.min(), WAVELENGTHS.max()+4.0
    with pytest.raises(RuntimeError):
        au._extrapolate_grid(WAVELENGTHS, wave_range, 1)


def test_extrapolate_bad_inputs():
    with pytest.raises(ValueError):
        au._extrapolate_grid(WAVELENGTHS, (2.9, 2.1))
    with pytest.raises(ValueError):
        au._extrapolate_grid(WAVELENGTHS, (4.1, 4.2))


def test_grid_from_map(wave_map, trace_profile):
    """Covers expected behavior of grid_from_map, including coverage of a previous bug
    where bad wavelengths were not being ignored properly"""

    wave_map = wave_map[0].copy()
    wavelengths = wave_map[0][::-1]
    trace_profile = trace_profile[0].copy()
    wave_grid = au.grid_from_map(wave_map, trace_profile, wave_range=None)

    assert np.allclose(wave_grid, wavelengths)

    # test custom wave_range
    wave_range = [wavelengths[2], wavelengths[-2]+0.01]
    wave_grid = au.grid_from_map(wave_map, trace_profile, wave_range=wave_range)
    assert np.allclose(wave_grid, wavelengths[2:-1])

    # test custom wave_range with extrapolation
    wave_range = [wavelengths[2], wavelengths[-1]+1]
    wave_grid = au.grid_from_map(wave_map, trace_profile, wave_range=wave_range)
    assert len(wave_grid) > len(wavelengths[2:])
    n_inside = wavelengths[2:].size
    assert np.allclose(wave_grid[:n_inside], wavelengths[2:])

    with pytest.raises(ValueError):
        au.grid_from_map(wave_map, trace_profile, wave_range=[0.1,0.2])


def xsinx(x):
    return x*np.sin(x)


def test_estim_integration_error():
    """
    Use as truth the x sin(x) from 0 to pi, has an analytic solution == pi.
    """

    n = 11
    grid = np.linspace(0, np.pi, n)
    err, rel_err = au._estim_integration_err(grid, xsinx)

    assert len(rel_err) == n-1
    assert np.all(rel_err >= 0)
    assert np.all(rel_err < 1)


@pytest.mark.parametrize("max_iter, rtol", [(1,1e-3), (10, 1e-9), (10, 1e-3), (1, 1e-9)])
def test_adapt_grid(max_iter, rtol):
    """
    Use as truth the x sin(x) from 0 to pi, has an analytic solution == pi.
    """

    input_grid = np.linspace(0, np.pi, 11)
    input_grid_diff = input_grid[1] - input_grid[0]
    max_grid_size = 100
    grid, is_converged = au._adapt_grid(input_grid,
                                        xsinx,
                                        max_grid_size,
                                        max_iter=max_iter,
                                        rtol=rtol)

    # ensure grid respects max_grid_size and max_iter in all cases
    assert len(grid) <= max_grid_size
    grid_diff = grid[1:] - grid[:-1]
    assert np.min(grid_diff) >= input_grid_diff/(2**max_iter)

    numerical_integral = np.trapz(xsinx(grid), grid)

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
        assert len(grid) <= 2*len(input_grid)


def test_adapt_grid_bad_inputs():
    with pytest.raises(ValueError):
        # input grid larger than max_grid_size
        au._adapt_grid(np.array([1,2,3]), xsinx, 2)


def test_trim_grids():

    grid_range = (-3, 3)
    grid0 = np.linspace(-3, 0, 4) # kept entirely.
    grid1 = np.linspace(-3, 0, 16) # removed entirely. Finer spacing doesn't matter, preceded by grid0
    grid2 = np.linspace(0, 3, 5) # kept from 0 to 3
    grid3 = np.linspace(-4, 4, 5) # removed entirely. Outside of grid_range and the rest is superseded

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
    grid0 = np.linspace(0, np.pi/2, 6) # kept entirely.
    grid1 = np.linspace(0, np.pi/2, 15) # removed entirely. Finer spacing doesn't matter, preceded by grid0
    grid2 = np.linspace(np.pi/2, np.pi, 11) # kept from pi/2 to pi
    
    # purposely make same lower index for grid2 as upper index for grid0 to test uniqueness of output

    all_grids = [grid0, grid1, grid2]
    all_estimate = [xsinx, xsinx, xsinx]

    rtol = 1e-3
    combined_grid = au.make_combined_adaptive_grid(all_grids, all_estimate, grid_range,
                                max_iter=10, rtol=rtol, max_total_size=100)
    
    numerical_integral = np.trapz(xsinx(combined_grid), combined_grid)

    assert np.unique(combined_grid).size == combined_grid.size
    assert np.isclose(numerical_integral, np.pi, rtol=rtol)


def test_throughput_soss():

    wavelengths = np.linspace(2,5,10)
    throughputs = np.ones_like(wavelengths)
    interpolator = au.ThroughputSOSS(wavelengths, throughputs)

    # test that it returns 1 for all wavelengths inside range
    interp = interpolator(wavelengths)
    assert np.allclose(interp[1:-1], throughputs[1:-1])
    assert interp[0] == 0
    assert interp[-1] == 0

    # test that it returns 0 for all wavelengths outside range
    wavelengths_outside = np.linspace(1,1.5,5)
    interp = interpolator(wavelengths_outside)
    assert np.all(interp == 0)

    # test ValueError raise for shape mismatch
    with pytest.raises(ValueError):
        au.ThroughputSOSS(wavelengths, throughputs[:-1])


def test_webb_kernel(webb_kernels, wave_map):

    wave_trace = wave_map[0][0]
    min_trace, max_trace = np.min(wave_trace), np.max(wave_trace)
    kern = webb_kernels[0]

    # basic ensure that the input is stored and shapes
    assert kern.wave_kernels.shape == kern.kernels.shape

    # test that pixels is mirrored around the center and has zero at center
    assert np.allclose(kern.pixels + kern.pixels[::-1], 0)
    assert kern.pixels[kern.pixels.size//2] == 0

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
    assert np.allclose(kern.poly[:,0], 0, atol=1e-1)
    assert np.allclose(kern.poly[:,1], kern.wave_center, atol=1e-1)

    # test interpolation function, which takes in a pixel and a wavelength and returns a throughput
    # this should return the triangle function at all wavelengths and zero outside range
    pix_half = kern.n_pix//2
    wl_test = np.linspace(min_trace, max_trace, 10)
    pixels_test = np.array([-pix_half-1, 0, pix_half, pix_half+1])
    
    data_in = kern.kernels[:,0]
    m = kern.min_value
    expected = np.array([m, np.max(data_in), m, m])

    interp = kern.f_ker(pixels_test, wl_test)
    assert interp.shape == (pixels_test.size, wl_test.size)
    diff = interp[:,1:] - interp[:,:-1]
    assert np.allclose(diff, 0)
    assert np.allclose(interp[:,0], expected, rtol=1e-3)

    # call the kernel object directly
    # this takes a wavelength and a central wavelength of the kernel, 
    # then converts to pixels to use self.f_ker internally
    assert kern(wl_test, wl_test).ndim == 1
    assert np.allclose(kern(wl_test, wl_test), np.max(data_in))

    # both inputs need to be same shape
    with pytest.raises(ValueError):
        kern(wl_test, wl_test[:-1])


def test_get_c_matrix(kernel):

    #TODO: what to use for grid? is it / can it be the same as wave_trace?
    # I think it can be the same but does not need to be, and would be a better
    # test if it were different, because the kernel is an interpolator that was
    # created using wavelengths that are included in wave_trace.
    matrix = au.get_c_matrix(kernel, grid, i_bounds=None, thresh=1e-5)

    # test with WebbKernel as the kernel
    # ensure normalized
    # ensure sparse


    # test where input kernel is a 2-D array instead of callable


    # test where input kernel is size 1


    # test where i_bounds is not None


    # Test invalid kernel input (wrong dimensions)


def test_finite_first_diff():

    wave_grid = np.linspace(0, 2*np.pi, 100)
    test_0 = np.ones_like(wave_grid)
    test_sin = np.sin(wave_grid)

    first_d = au.finite_first_d(wave_grid)
    assert first_d.size == (wave_grid.size - 1)*2

    # test trivial example returning zeros for constant
    f0 = first_d.dot(test_0)
    assert np.allclose(f0, 0)

    # test derivative of sin returns cos
    wave_between = (wave_grid[1:] + wave_grid[:-1])/2
    f_sin = first_d.dot(test_sin)
    assert np.allclose(f_sin, np.cos(wave_between), atol=1e-3)


def test_extraction_engine_init(
    wave_map,
    trace_profile,
    throughput,
    wave_grid,
    mask_trace_profile,
    detector_mask,
    engine,
):
    """Test the init of the engine with default/good inputs"""
    
    # test wave_grid became unique
    assert engine.wave_grid.dtype == np.float64
    unq = np.unique(wave_grid)
    assert np.allclose(engine.wave_grid, unq)
    assert engine.n_wavepoints == unq.size

    for order in [0,1]:
        # test assignment of attributes and conversion to expected float64 dtype
        assert engine.wave_map[order].dtype == np.float64
        assert engine.trace_profile[order].dtype == np.float64
        assert engine.kernels[order].dtype == np.float64
        assert engine.mask_trace_profile[order].dtype == np.bool_
        
        assert np.allclose(engine.wave_map[order], wave_map[order])
        assert np.allclose(engine.trace_profile[order], trace_profile[order])
        assert np.allclose(engine.mask_trace_profile[order], mask_trace_profile[order])

        # test derived attributes
        assert engine.data_shape == DATA_SHAPE
        assert engine.n_orders == 2
        
    # test wave_p and wave_m. separate unit test for their calculation
    for att in ["wave_p", "wave_m"]:
        wave = getattr(engine, att)
        assert wave.dtype == np.float64
        assert wave.shape == (2,)+DATA_SHAPE

    # test _get_i_bounds
    assert len(engine.i_bounds) == 2
    for order in [0,1]:
        assert len(engine.i_bounds[order]) == 2
        assert engine.i_bounds[order][0] >= 0
        assert engine.i_bounds[order][1] < DATA_SHAPE[1]
        assert all([isinstance(val, int) for val in engine.i_bounds[order]])

    # test to ensure that wave_map is considered for bounds
    # in order 1 the wave_map is more restrictive on the shortwave end
    # in order 2 the wave_grid is more restrictive on the longwave end
    # in order 1 no restriction on the longwave end so get the full extent
    # in order 2 no restriction on the shortwave end so get the full extent
    assert engine.i_bounds[0][0] > 0
    assert engine.i_bounds[1][1] < engine.n_wavepoints
    # TODO: off-by-one error here. why does this fail?
    # check what this looks like on a real run on main
    # assert engine.i_bounds[0][1] == engine.n_wavepoints
    # assert engine.i_bounds[1][0] == 0


    # test _get_masks
    # ensure they all include the input detector bad pixel mask
    for mask in [engine.mask_ord[0], engine.mask_ord[1], engine.mask, engine.general_mask]:
        assert mask.dtype == np.bool_
        assert mask.shape == DATA_SHAPE
        assert np.all(mask[detector_mask == 1])

    # general_mask should be a copy of mask
    assert np.allclose(engine.mask, engine.general_mask)

    wave_bnds = [WAVE_BNDS_O1, WAVE_BNDS_O2]
    for order in [0,1]:
        # ensure wavelength bounds from wave_grid are respected in mask_ord
        mask = engine.mask_ord[order]
        wls = wave_map[order]
        lo, hi = wave_bnds[order]
        outside = (wls > lo) | (wls < hi)
        assert np.all(mask[outside])

        # a bit paradoxically, engine.mask_ord does not contain a single order's mask_trace_profile
        # instead, it's mask_trace_profile[0] AND mask_trace_profile[1], i.e., 
        # the trace profiles of BOTH orders are UNmasked in the masks of each order
        # the general mask and wavelength bounds are then applied, so the
        # only difference between mask_ord[0] and mask_ord[1] are the wavelength bounds
        # test that NOT all locations masked by mask_trace_profile are masked in mask_ord
        assert not np.all(mask[mask_trace_profile[order]])
        # test that all locations masked by both profiles are masked in mask_ord
        combined_profile = mask_trace_profile[0] & mask_trace_profile[1]
        assert np.all(mask[combined_profile])

    # test throughput function conversion to array
    for order in [0,1]:
        thru = engine.throughput[order]
        assert thru.size == engine.n_wavepoints
        assert np.all(thru >= 0)
        assert thru.dtype == np.float64

    # test kernel is cast to proper shape for input (trivial) kernel
    for order in [0,1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        expected_shape = (n_valid, engine.n_wavepoints)
        assert engine.kernels[order].shape == expected_shape
        # for trivial kernel only one element per row is nonzero
        assert engine.kernels[order].count_nonzero() == expected_shape[0]

    # test weights. see separate unit tests to ensure the calculation is correct
    for order in [0,1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        weights = engine.weights[order]
        k_idx = engine.weights_k_idx[order]
        assert weights.dtype == np.float64
        assert np.issubdtype(k_idx.dtype, np.integer)

        # TODO: why is weights the same size as ~engine.mask, and not ~engine.mask_ord?
        assert weights.shape == (np.sum(~engine.mask), n_valid)
        assert k_idx.shape[0] == weights.shape[0]

    # test assignment of empty attributes
    for att in ["w_t_wave_c", "tikho_mat", "_tikho_mat"]:
        assert hasattr(engine, att)
    assert getattr(engine, "pixel_mapping") == [None, None]


def test_extraction_engine_bad_inputs(
    wave_map,
    trace_profile,
    throughput,
    kernels_unity,
    wave_grid,
    mask_trace_profile,
    detector_mask,
):
    # not enough good pixels in order
    with pytest.raises(atoca.MaskOverlapError):
        detector_mask = np.ones_like(detector_mask)
        detector_mask[5:7,50:55] = 0 #still a few good pixels but very few
        atoca.ExtractionEngine(wave_map,
                               trace_profile,
                               throughput,
                               kernels_unity,
                               wave_grid,
                               mask_trace_profile,
                               global_mask=detector_mask)

    # wrong number of orders
    with pytest.raises(ValueError):
        atoca.ExtractionEngine(wave_map,
                               trace_profile,
                               throughput,
                               kernels_unity,
                               wave_grid,
                               mask_trace_profile,
                               global_mask=detector_mask,
                               orders=[0,])


def test_get_attributes(engine):
    # test string input
    assert np.allclose(engine.get_attributes("wave_map"), engine.wave_map)

    # test list of strings input
    name_list = ["wave_map", "wave_grid"]
    att_list = engine.get_attributes(*name_list)
    expected = [engine.wave_map, engine.wave_grid]
    for i in range(len(expected)):
        for j in range(2): #orders
            assert np.allclose(att_list[i][j], expected[i][j])

    # test i_order not None
    att_list = engine.get_attributes(*name_list, i_order=1) 
    expected = [engine.wave_map[1], engine.wave_grid[1]]
    for i in range(len(expected)):
        assert np.allclose(att_list[i], expected[i])


def test_update_throughput(engine, throughput):

    old_thru = engine.throughput

    # test callable input
    new_thru = [throughput[1], throughput[1]]
    engine.update_throughput(new_thru)
    for i, thru in enumerate(engine.throughput):
        assert isinstance(thru, np.ndarray)
        assert thru.shape == engine.wave_grid.shape
    assert np.allclose(engine.throughput[0], engine.throughput[1])

    # test array input
    # first reset to old throughput
    engine.throughput = old_thru
    new_thru = [thru*2 for thru in engine.throughput]
    engine.update_throughput(new_thru)
    for i, thru in enumerate(engine.throughput):
        assert np.allclose(thru, old_thru[i]*2)

    # test fail on bad array shape
    new_thru = [thru[:-1] for thru in engine.throughput]
    with pytest.raises(ValueError):
        engine.update_throughput(new_thru)

    # test fail on callable that doesn't return correct array shape
    def new_thru_f(wl):
        return 1.0
    with pytest.raises(ValueError):
        engine.update_throughput([new_thru_f, new_thru_f])


def test_create_kernels(engine):
    # TODO: make a non-trivial kernel fixture to work with this
    # trivial example already done
    pass


def test_wave_grid_c(engine):
    for order in [0,1]:
        n_valid = engine.i_bounds[order][1] - engine.i_bounds[order][0]
        assert engine.wave_grid_c(order).size == n_valid


def test_set_w_t_wave_c(engine):
    """all this does is copy whatever is input"""
    product = np.zeros((1,))
    engine._set_w_t_wave_c(0, product)
    assert len(engine.w_t_wave_c) == engine.n_orders
    assert engine.w_t_wave_c[0] == product
    assert engine.w_t_wave_c[1] == []
    assert product is not engine.w_t_wave_c[0]


def test_get_pixel_mapping(engine):

    pixel_mapping_0 = engine.get_pixel_mapping(0)
    # check attribute is set and identical to output
    # check the second one is not set but there is space for it
    assert hasattr(engine, "pixel_mapping")
    assert len(engine.pixel_mapping) == engine.n_orders
    assert np.allclose(engine.pixel_mapping[0].data, pixel_mapping_0.data)
    assert engine.pixel_mapping[1] is None

    # set the second one so can check both at once
    engine.get_pixel_mapping(1)
    for order in [0,1]:
        mapping = engine.pixel_mapping[order]
        assert mapping.dtype == np.float64
        # TODO: why is this the shape, instead of using mask_ord and only the valid wave_grid?
        expected_shape = (np.sum(~engine.mask), engine.wave_grid.size)
        assert mapping.shape == expected_shape

        # test that w_t_wave_c is getting saved
        w_t_wave_c = engine.w_t_wave_c[order]
        assert w_t_wave_c.dtype == np.float64
        assert w_t_wave_c.shape == expected_shape
    
        # check if quick=True works
        mapping_quick = engine.get_pixel_mapping(order, quick=True)
        assert np.allclose(mapping.data, mapping_quick.data)

    # check that quick=True does not work if w_t_wave_c unsaved
    engine.w_t_wave_c = None
    with pytest.raises(AttributeError):
        engine.get_pixel_mapping(1, quick=True)

    # TODO: follow the shape of every sparse matrix through this function, see if
    # it makes sense or if it would make more sense to mask differently
    # TODO: test that the math is actually correct using a trivial example (somehow)


def test_rebuild(engine):

    detector_model = engine.rebuild(f_lam)
    assert detector_model.dtype == np.float64
    assert detector_model.shape == engine.wave_map[0].shape

    # test that input spectrum is ok as either callable or array
    assert np.allclose(detector_model, engine.rebuild(f_lam(engine.wave_grid)))

    # test fill value
    detector_model_nans = engine.rebuild(f_lam, fill_value=np.nan)
    assert np.allclose(np.isnan(detector_model_nans), engine.general_mask)


def f_lam(wl, m=SPECTRAL_SLOPE, b=0):
    """
    Estimator for flux as function of wavelength
    Returns linear function of wl with slope m and intercept b
    
    This function is also used in this test suite as 
    """
    return m*wl + b


@pytest.fixture(scope="module")
def imagemodel(engine, detector_mask):
    """
    use engine.rebuild to make an image model from an expected f(lambda).
    Then we can ensure it round-trips
    """

    rng = np.random.default_rng(seed=42)
    shp = engine.trace_profile[0].shape

    # make the detector bad values NaN, but leave the trace masks alone
    # in reality, of course, this is backward: detector bad values
    # would be determined from data
    data = engine.rebuild(f_lam, fill_value=0.0)
    data[detector_mask] = np.nan

    # add noise
    noise_scaling = 3e-5
    data += noise_scaling*rng.standard_normal(shp)

    # error random, all positive, but also always larger than a certain value
    # to avoid very large values of data / error
    error = noise_scaling*(rng.standard_normal(shp)**2 + 0.5)

    return data, error


def test_build_sys(imagemodel, engine):
    
    data, error = imagemodel
    matrix, result = engine.build_sys(data, error)
    assert result.size == engine.n_wavepoints
    assert matrix.shape == (result.size, result.size)



def test_get_detector_model(imagemodel, engine):

    data, error = imagemodel
    unmasked_size = np.sum(~engine.mask)
    b_matrix, data_matrix = engine.get_detector_model(data, error)

    assert data_matrix.shape == (1, unmasked_size)
    assert b_matrix.shape == (unmasked_size, engine.n_wavepoints)
    assert np.allclose(data_matrix.toarray()[0], (data/error)[~engine.mask])


def test_estimate_tikho_factors(engine):
    
    factor = engine.estimate_tikho_factors(f_lam)
    assert isinstance(factor, float)


    # very approximate calculation of tik fac looks like
    # n_pixels = (~engine.mask).sum()
    # flux = f_lam(engine.wave_grid)
    # dlam = engine.wave_grid[1:] - engine.wave_grid[:-1]
    # print(n_pixels/np.mean(flux[1:] * dlam))


@pytest.fixture(scope="module")
def tikho_tests(imagemodel, engine):
    data, error = imagemodel

    log_guess = np.log10(engine.estimate_tikho_factors(f_lam))
    factors = np.logspace(log_guess - 9, log_guess + 9, 19)
    return factors, engine.get_tikho_tests(factors, data, error)


def test_get_tikho_tests(tikho_tests, engine):

    factors, tests = tikho_tests
    unmasked_size = np.sum(~engine.mask)

    # test all the output shapes
    assert np.allclose(tests["factors"], factors)
    assert tests["solution"].shape == (len(factors), engine.n_wavepoints)
    assert tests["error"].shape == (len(factors), unmasked_size)
    assert tests["reg"].shape == (len(factors), engine.n_wavepoints-1)
    assert tests["chi2"].shape == (len(factors),)
    assert tests["chi2_soft_l1"].shape == (len(factors),)
    assert tests["chi2_cauchy"].shape == (len(factors),)
    assert np.allclose(tests["grid"], engine.wave_grid)

    # test data type is preserved through solve
    for key in tests.keys():
        assert tests[key].dtype == np.dtype("float64")


def test_best_tikho_factor(engine, tikho_tests):

    input_factors, tests = tikho_tests
    fit_modes = ["all", "curvature", "chi2", "d_chi2"]
    best_factors = []
    for mode in fit_modes:
        factor = engine.best_tikho_factor(tests, mode)
        assert isinstance(factor, float)
        best_factors.append(factor)
    
    # ensure fit_mode=all found one of the three others
    assert best_factors[0] in best_factors[1:]

    # TODO: test the logic tree by manually changing the tests dict
    # this is non-trivial because dchi2 and curvature
    # are both derived from other keys in the dictionary, not just the statistical metrics


def test_call(engine, tikho_tests, imagemodel):
    """
    Run the actual extract method.
    Ensure it can retrieve the input spectrum based on f_lam to within a few percent
    at all points on the wave_grid.

    Note this round-trip implicitly checks the math of the build_sys, get_detector_model,
    _solve, and _solve_tikho, at least at first blush.
    """
    data, error = imagemodel
    _, tests = tikho_tests
    best_factor = engine.best_tikho_factor(tests, "all")

    expected_spectrum = f_lam(engine.wave_grid)
    for tikhonov in [True, False]:
        spectrum = engine(data, error, tikhonov=tikhonov, factor=best_factor)
        diff = (spectrum - expected_spectrum)/expected_spectrum
        assert not np.all(np.isnan(diff))
        diff = diff[~np.isnan(diff)]
        assert np.all(np.abs(diff) < 0.05)
    
    # test bad input, failing to put factor in for Tikhonov solver
    with pytest.raises(ValueError):
        engine(data, error, tikhonov=True)


def test_compute_likelihood(engine, imagemodel):
    """Ensure log-likelihood is highest for the correct slope"""

    data, error = imagemodel
    test_slopes = np.arange(0, 5, 0.5)
    logl = []
    for slope in test_slopes:
        spectrum = partial(f_lam, m=slope)
        logl.append(engine.compute_likelihood(spectrum, data, error))
    
    assert np.argmax(logl) == np.argwhere(test_slopes == SPECTRAL_SLOPE)
    assert np.all(np.array(logl) < 0)


def test_sparse_c(kern_array):
    """Here kernel must be a 2-D array already, of shape (N_ker, N_k_convolved)"""

    # test typical case n_k = n_kc and i=0
    n_k = kern_array.shape[1]
    i_zero = 0
    matrix = au._sparse_c(kern_array, n_k, i_zero)

    # TODO: add more here


@pytest.fixture(scope="module")
def tikhoTests():
    """Make a TikhoTests dictionary"""


    return au.TikhoTests({'factors': factors,
                          'solution': sln,
                          'error': err,
                          'reg': reg,
                          'grid': wave_grid})



def test_tikho_tests(tikhoTests):

    assert False