import numpy as np
import pytest

from jwst.extract_1d.soss_extract.pastasoss import (
    _extrapolate_to_wavegrid,
    _find_spectral_order_index,
    _get_soss_traces,
    _get_wavelengths,
    get_soss_traces,
    get_soss_wavemaps,
)
from jwst.extract_1d.soss_extract.tests.helpers import (
    PWCPOS,
    TRACE_END_IDX,
    WAVE_BNDS_O1,
    WAVE_BNDS_O2,
)

"""Test coverage for the helper functions in pastasoss.py"""


def test_wavecal_models(refmodel):
    wave_bnds = [WAVE_BNDS_O1, WAVE_BNDS_O2]
    for order in [1, 2]:
        idx = order - 1
        bnds = wave_bnds[idx]
        x = np.arange(0, TRACE_END_IDX[idx] + 1)
        wavelengths = _get_wavelengths(refmodel, x, PWCPOS, order)

        # check shapes
        assert wavelengths.shape == x.shape
        assert np.isclose(wavelengths[0], bnds[0])
        assert np.isclose(wavelengths[-1], bnds[1])

        # ensure unique and descending
        diff = wavelengths[1:] - wavelengths[:-1]
        assert np.all(diff < 0)


def test_rotate():
    # TODO: add meaningful tests of rotate
    pass


def test_find_spectral_order_index(refmodel):
    """TODO: why doesn't this raise an error when order is not recognized?
    Surely it's a bad idea to have the index set to -1?"""
    for order in [1, 2]:
        idx = _find_spectral_order_index(refmodel, order)
        assert idx == order - 1

    for order in [0, "bad", None]:
        with pytest.raises(ValueError):
            _find_spectral_order_index(refmodel, order)


def test_get_soss_traces(refmodel):
    for order in ["1", "2"]:
        idx = int(order) - 1
        for subarray in ["SUBSTRIP96", "SUBSTRIP256"]:
            order_out, x_new, y_new, wavelengths = _get_soss_traces(
                refmodel, PWCPOS, order, subarray
            )

            assert str(order_out) == order
            # since always interpolated back to original x, x_new should equal x
            x_in, y_in = refmodel.traces[idx].trace.T.copy()
            assert np.allclose(x_new, x_in)
            # and wavelengths are same as what you get from _get_wavelengths on x_in
            wave_expected = _get_wavelengths(refmodel, x_in, PWCPOS, int(order))
            assert np.allclose(wavelengths, wave_expected)

            # the y coordinate is the tricky one. it was rotated by pwcpos - refmodel.meta.pwcpos_cmd
            # about pivot_x, pivot_y
            assert y_new.shape == wavelengths.shape
            # TODO: add meaningful tests of y


def test_extrapolate_to_wavegrid(refmodel):
    wavemin = 0.5
    wavemax = 5.5
    nwave = 501
    wave_grid = np.linspace(wavemin, wavemax, nwave)

    # only test first order
    x = np.arange(0, TRACE_END_IDX[0] + 1)
    wl = _get_wavelengths(refmodel, x, PWCPOS, 1)

    # first ensure test setup gives all wl in wave_grid
    # floating-point precision issues make np.around calls necessary
    assert np.all(np.isin(np.around(wl, 5), np.around(wave_grid, 5)))

    x_extrap = _extrapolate_to_wavegrid(wave_grid, wl, x)
    assert x_extrap.shape == wave_grid.shape

    # test that all x in x_extrap
    assert np.all(np.isin(np.around(x, 5), np.around(x_extrap, 5)))

    # test extrapolated slope is same as input slope, since these are linear
    m_extrap = (x_extrap[-1] - x_extrap[0]) / (wave_grid[-1] - wave_grid[0])
    m = (x[-1] - x[0]) / (wl[-1] - wl[0])
    assert np.isclose(m_extrap, m)


@pytest.mark.parametrize("pwcpos", [245.79 - 0.24, 245.79, 245.79 + 0.24])  # edges of bounds
@pytest.mark.parametrize("order", ["1", "2"])
@pytest.mark.parametrize("subarray", [None, "SUBSTRIP96", "FULL"])
def test_get_soss_traces_public(subarray, order, pwcpos):
    """Test of public interface to get_soss_traces, which should not require datamodel or refmodel"""
    order_out, x_new, y_new, wavelengths = get_soss_traces(pwcpos, order, subarray=subarray)

    assert str(order_out) == order
    assert x_new.shape == y_new.shape
    assert x_new.shape == wavelengths.shape


@pytest.mark.parametrize("pwcpos", [245.79 - 0.24, 245.79, 245.79 + 0.24])  # edges of bounds
@pytest.mark.parametrize("padsize", [None, 9])
@pytest.mark.parametrize("subarray", ["SUBSTRIP256", "SUBSTRIP96", "FULL"])
def test_get_soss_wavemaps_public(subarray, padsize, pwcpos):
    """Test of public interface to get_soss_wavemaps, which should not require datamodel or refmodel"""
    subarray_shapes = {"SUBSTRIP96": 96, "SUBSTRIP256": 256, "FULL": 2048}
    wavemaps, traces = get_soss_wavemaps(
        pwcpos, subarray=subarray, padsize=padsize, spectraces=True
    )

    if padsize is None:
        padsize = 0
    expected_shape = (2, subarray_shapes[subarray] + padsize * 2, 2048 + padsize * 2)
    assert wavemaps.shape == expected_shape
    assert traces.shape == (2, 3, 5001)


def test_get_soss_traces_bad_pwcpos():
    """Test that get_soss_traces raises ValueError for bad PWC positions"""
    with pytest.raises(ValueError, match="PWC position 245.0 is outside bounds"):
        get_soss_traces(245.0, "1")
    with pytest.raises(ValueError, match="PWC position 247.0 is outside bounds"):
        get_soss_traces(247.0, "2")
