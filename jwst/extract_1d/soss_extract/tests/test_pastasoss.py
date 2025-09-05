import numpy as np
import pytest

from jwst.extract_1d.soss_extract.pastasoss import (
    _convert_refmodel_poly_to_astropy,
    _extrapolate_to_wavegrid,
    _find_spectral_order_index,
    _get_soss_traces,
    _get_wavelengths,
    _verify_requested_orders,
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


def test_verify_requested_orders(log_watcher):
    refmodel_orders = [1, 2, 3]
    requested_orders = [1, 2]
    good_orders = _verify_requested_orders(refmodel_orders, requested_orders)
    assert good_orders == requested_orders

    requested_orders = [1, 2, 4]
    watcher = log_watcher(
        "jwst.extract_1d.soss_extract.pastasoss",
        message="Some requested orders were not found in reference model",
    )
    good_orders = _verify_requested_orders(refmodel_orders, requested_orders)
    watcher.assert_seen()
    assert good_orders == [1, 2]

    with pytest.raises(ValueError):
        requested_orders = [0]
        _verify_requested_orders(refmodel_orders, requested_orders)


@pytest.mark.parametrize("degree", [3, 5])
def test_convert_refmodel_poly_to_astropy(degree):
    """Test that the reference model polynomial coefficients are correctly converted to Astropy format."""

    n_coeffs = (degree + 1) * (degree + 2) // 2
    coeffs = np.arange(n_coeffs)
    astropy_poly = _convert_refmodel_poly_to_astropy(coeffs)

    # pick some random numbers
    x = np.linspace(0, 1, 25)
    offset = np.linspace(-3, 3, 25)

    # this is what used to be in the code for degree 5
    poly_features = np.array(
        [
            x,
            offset,
            x**2,
            x * offset,
            offset**2,
            x**3,
            x**2 * offset,
            x * offset**2,
            offset**3,
            x**4,
            x**3 * offset,
            x**2 * offset**2,
            x * offset**3,
            offset**4,
            x**5,
            x**4 * offset,
            x**3 * offset**2,
            x**2 * offset**3,
            x * offset**4,
            offset**5,
        ]
    )

    poly_features = poly_features[: n_coeffs - 1, :]
    expected = coeffs[0] + coeffs[1:] @ poly_features
    astropy_result = astropy_poly(x, offset)
    assert np.allclose(astropy_result, expected)


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
    for order in ["1", "2", "3"]:
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
@pytest.mark.parametrize("orders_requested", [[1], [2], [1, 2]])
def test_get_soss_wavemaps_public(subarray, padsize, pwcpos, orders_requested):
    """Test of public interface to get_soss_wavemaps, which should not require datamodel or refmodel"""
    subarray_shapes = {"SUBSTRIP96": 96, "SUBSTRIP256": 256, "FULL": 2048}
    wavemaps, traces = get_soss_wavemaps(
        pwcpos,
        subarray=subarray,
        padsize=padsize,
        spectraces=True,
        orders_requested=orders_requested,
    )
    n_orders = len(orders_requested)
    if padsize is None:
        padsize = 0
    expected_shape = (n_orders, subarray_shapes[subarray] + padsize * 2, 2048 + padsize * 2)
    assert wavemaps.shape == expected_shape
    assert traces.shape == (n_orders, 3, 5001)


def test_get_soss_traces_bad_pwcpos():
    """Test that get_soss_traces raises ValueError for bad PWC positions"""
    with pytest.raises(ValueError, match="PWC position 245.0 is outside bounds"):
        get_soss_traces(245.0, "1")
    with pytest.raises(ValueError, match="PWC position 247.0 is outside bounds"):
        get_soss_traces(247.0, "2")
