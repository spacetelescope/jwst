
import pytest
from jwst.extract_1d.soss_extract import atoca_utils as au
import numpy as np

def test_arange_2d():

    starts = np.array([3,4,5])
    stops = np.ones(starts.shape)*7
    out = au.arange_2d(starts, stops)

    bad = 65535
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


# wavelengths have min, max (1.5, 4.0) and a bit of non-linearity
WAVELENGTHS = np.linspace(1.5, 3.0, 50) + np.sin(np.linspace(0, np.pi/2, 50))
@pytest.fixture(scope="function")
def wave_map():
    wave_map = np.array([
        WAVELENGTHS,
        WAVELENGTHS+0.2,
        WAVELENGTHS+0.2,
        WAVELENGTHS+0.2,
        WAVELENGTHS+0.4,
    ])
    wave_map[1,3] = -1 #test skip of bad value
    return wave_map


@pytest.fixture(scope="function")
def wave_map_o2(wave_map):
    return np.copy(wave_map) - 1.0
    

@pytest.fixture(scope="function")
def trace_profile(wave_map):
    thrpt = np.array([0.01, 0.95, 1.0, 0.8, 0.01])
    trace_profile = np.ones_like(wave_map)
    return trace_profile*thrpt[:,None]


@pytest.fixture(scope="function")
def trace_profile_o2(wave_map_o2):
    thrpt = np.array([0.001, 0.01, 0.01, 0.2, 0.99])
    trace_profile = np.ones_like(wave_map_o2)
    return trace_profile*thrpt[:,None]


@pytest.mark.parametrize("dispersion_axis", [0,1])
def test_get_wv_map_bounds(wave_map, dispersion_axis):
    """
    top is the low-wavelength end, bottom is high-wavelength end
    """
    if dispersion_axis == 0:
        wave_flip = wave_map.T
    else:
        wave_flip = wave_map
    wave_top, wave_bottom = au._get_wv_map_bounds(wave_flip, dispersion_axis=dispersion_axis)

    # flip the results back so we can re-use the same tests
    if dispersion_axis == 0:
        wave_top = wave_top.T
        wave_bottom = wave_bottom.T

    diff = (WAVELENGTHS[1:]-WAVELENGTHS[:-1])/2
    diff_lower = np.insert(diff,0,diff[0])
    diff_upper = np.append(diff,diff[-1])
    wave_top_expected = WAVELENGTHS-diff_lower
    wave_bottom_expected = WAVELENGTHS+diff_upper

    # basic test
    assert wave_top.shape == wave_bottom.shape == (wave_map.shape[0],)+WAVELENGTHS.shape
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


@pytest.mark.parametrize("wave_range", [(2.1, 3.9), (1.8, 4.5)])
def test__extrapolate_grid(wave_range):

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

    wave_grid = au.grid_from_map(wave_map, trace_profile, wave_range=None)

    # expected output is very near WAVELENGTHS+0.2 because that's what all the high-weight
    # rows of the wave_map are set to.
    assert np.allclose(wave_grid, wave_map[2])

    # test custom wave_range
    wave_range = [wave_map[2,2], wave_map[2,-2]+0.01]
    wave_grid = au.grid_from_map(wave_map, trace_profile, wave_range=wave_range)
    assert np.allclose(wave_grid, wave_map[2,2:-1])

    # test custom wave_range with extrapolation
    wave_range = [wave_map[2,2], wave_map[2,-1]+1]
    wave_grid = au.grid_from_map(wave_map, trace_profile, wave_range=wave_range)
    assert len(wave_grid) > len(wave_map[2,2:])
    n_inside = wave_map[2,2:].size
    assert np.allclose(wave_grid[:n_inside], wave_map[2,2:])

    with pytest.raises(ValueError):
        au.grid_from_map(wave_map, trace_profile, wave_range=[0.5,0.9])


@pytest.mark.parametrize("n_os", ([4,1], 1))
def test__get_soss_grid(n_os, wave_map, trace_profile, wave_map_o2, trace_profile_o2):
    """
    wave_map has min, max wavelength of 1.5, 4.0 but throughput makes this 1.7, 4.2
    wave_map_o2 has min, max wavelength of 0.5, 3.0 but throughput makes this 0.9, 3.4
    """

    # choices of wave_min, wave_max force extrapolation on low end 
    # and also makes sure both orders matter
    wave_min = 0.55
    wave_max = 4.0
    wave_maps = np.array([wave_map, wave_map_o2])
    trace_profiles = np.array([trace_profile, trace_profile_o2])
    wave_grid = au._get_soss_grid(wave_maps, trace_profiles, wave_min, wave_max, n_os)

    delta_lower = wave_grid[1]-wave_grid[0]
    delta_upper = wave_grid[-1]-wave_grid[-2]

    # ensure no duplicates and strictly ascending
    assert wave_grid.size == np.unique(wave_grid).size
    assert np.all(wave_grid[1:] > wave_grid[:-1])

    # esnure grid is within bounds but just abutting bounds
    assert wave_grid.min() >= wave_min
    assert wave_grid.min() <= wave_min + delta_lower
    assert wave_grid.max() <= wave_max
    assert wave_grid.max() >= wave_max - delta_upper

    # ensure oversample factor changes for different n_os
    # this is a bit complicated because the wavelength spacing is nonlinear in wave_map
    # by a factor of 2 to begin with, so check that the ratio of the wavelength spacing
    # in the upper vs lower end of the wl ranges look approx like what was input, modulo
    # the oversample factor of the two orders
    if n_os == 1:
        n_os = [1,1]
    og_spacing_lower = WAVELENGTHS[1]-WAVELENGTHS[0]
    og_spacing_upper = WAVELENGTHS[-1]-WAVELENGTHS[-2]
    expected_ratio = int((n_os[0]/n_os[1])*(np.around(og_spacing_lower/og_spacing_upper)))

    spacing_lower = np.mean(wave_grid[1:6]-wave_grid[:5])
    spacing_upper = np.mean(wave_grid[-5:]-wave_grid[-6:-1])
    actual_ratio = int(np.around((spacing_lower/spacing_upper)))
    
    # for n=1 case we expect 2, for n=(4,1) case we expect 8
    assert expected_ratio == actual_ratio


def test__get_soss_grid_bad_inputs(wave_map, trace_profile):
    with pytest.raises(ValueError):
        # test bad input shapes
        au._get_soss_grid(wave_map, trace_profile, 0.5, 0.9, 1)

    wave_maps = np.array([wave_map, wave_map])
    trace_profiles = np.array([trace_profile, trace_profile])

    with pytest.raises(ValueError):
        # test bad n_os shape
        au._get_soss_grid(wave_maps, trace_profiles, 0.5, 0.9, [1,1,1])



