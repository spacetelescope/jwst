
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


FIBONACCI = np.array([1,1,2,3,5,8,13], dtype=float)
@pytest.fixture(scope="module")
def wave_map():
    wave_map = np.array([
        FIBONACCI,
        FIBONACCI+1,

    ])
    wave_map[1,3] = -1 #test skip of bad value
    return wave_map


@pytest.mark.parametrize("dispersion_axis", [0,1])
def test_get_wv_map_bounds(wave_map, dispersion_axis):
    """
    top is the low-wavelength end, bottom is high-wavelength end
    """
    if dispersion_axis == 0:
        wave_flip = wave_map.T
    else:
        wave_flip = wave_map
    wave_top, wave_bottom = au.get_wv_map_bounds(wave_flip, dispersion_axis=dispersion_axis)

    # flip the results back so we can re-use the same tests
    if dispersion_axis == 0:
        wave_top = wave_top.T
        wave_bottom = wave_bottom.T

    diff = (FIBONACCI[1:]-FIBONACCI[:-1])/2
    diff_lower = np.insert(diff,0,diff[0])
    diff_upper = np.append(diff,diff[-1])
    wave_top_expected = FIBONACCI-diff_lower
    wave_bottom_expected = FIBONACCI+diff_upper

    # basic test
    assert wave_top.shape == wave_bottom.shape == (2,)+FIBONACCI.shape
    assert np.allclose(wave_top[0], wave_top_expected)
    assert np.allclose(wave_bottom[0], wave_bottom_expected)

    # test skip bad pixel
    assert wave_top[1,3] == 0
    assert wave_bottom[1,3] == 0

    # test bad input error raises
    with pytest.raises(ValueError):
        au.get_wv_map_bounds(wave_flip, dispersion_axis=2)



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


