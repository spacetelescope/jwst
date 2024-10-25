
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
def test_get_wave_p_or_m(dispersion_axis):
    return