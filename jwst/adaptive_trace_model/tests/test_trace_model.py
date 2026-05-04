import numpy as np
import pytest
from stdatamodels.jwst.datamodels import ImageModel

from jwst.adaptive_trace_model import trace_model as tm
from jwst.adaptive_trace_model.tests import helpers


@pytest.fixture(scope="module")
def nrs_slit_model():
    model = helpers.nirspec_slit_model_with_source()
    yield model
    model.close()


@pytest.fixture(scope="module")
def fit_2d_spline_input(nrs_slit_model):
    slit = nrs_slit_model.slits[0]
    flux = slit.data
    wcs = slit.meta.wcs
    x, y = np.meshgrid(np.arange(flux.shape[1]), np.arange(flux.shape[0]))
    _, alpha, _ = wcs.transform("detector", "slit_frame", x, y)
    return slit, flux, alpha


def test_fit_2d_spline_trace(fit_2d_spline_input):
    slit, flux, alpha = fit_2d_spline_input
    splines = tm.fit_2d_spline_trace(flux, alpha)

    assert len(splines) == flux.shape[1]
    scales = [s["scale"] for s in splines.values()]

    # for unscaled fits, output scales should be close to 1
    np.testing.assert_allclose(list(scales), 1.0, atol=0.06)

    # evaluated fits should be close to input data
    region_map = (~np.isnan(slit.wavelength)).astype(int)
    trace_used, full_trace = tm._trace_image(flux.shape, {1: splines}, region_map, alpha)

    # Fit values should be close to flux
    atol = 0.25 * np.nanmax(flux)
    full_indx = (region_map == 1) & ~np.isnan(full_trace)
    assert np.sum(full_indx) > 0
    np.testing.assert_allclose(flux[full_indx], full_trace[full_indx], atol=atol)

    # Trace used is a subset of the full trace
    used_indx = (region_map == 1) & ~np.isnan(trace_used)
    assert 0 < np.sum(used_indx) < np.sum(full_indx)
    np.testing.assert_allclose(flux[used_indx], trace_used[used_indx], atol=atol)


def test_fit_2d_spline_trace_fail(monkeypatch, caplog, fit_2d_spline_input):
    # mock a failure in the fit
    def mock_fit(*args, **kwargs):
        raise RuntimeError("Bad fit")

    monkeypatch.setattr(tm, "bspline_fit", mock_fit)

    slit, flux, alpha = fit_2d_spline_input
    splines = tm.fit_2d_spline_trace(flux, alpha)
    assert len(splines) == 0
    assert "Spline fit failed" in caplog.text


def test_fit_2d_spline_trace_none(monkeypatch, fit_2d_spline_input):
    """Test that fit failure triggers a single retry with fewer knots"""

    # mock a quiet failure in the fit: always return None
    class MockFit(object):
        def __init__(self):
            self.call_count = 0

        def __call__(self, *args, **kwargs):
            # Track how many times this mock was called
            self.call_count += 1
            return None

    mock_fit = MockFit()
    monkeypatch.setattr(tm, "bspline_fit", mock_fit)

    slit, flux, alpha = fit_2d_spline_input
    splines = tm.fit_2d_spline_trace(flux, alpha)
    assert len(splines) == 0

    # The fit is re-attempted once at every column
    assert mock_fit.call_count == 2 * flux.shape[1]


@pytest.mark.parametrize(
    "mode, detector, slit, expected",
    [
        (
            "NRS_IFU",
            "NRS1",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "require_ngood": 15,
                "spline_bkpt": 68,
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "NRS_IFU",
            "NRS2",
            None,
            {
                "lrange": 50,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "require_ngood": 15,
                "spline_bkpt": 68,
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "NRS_SLIT",
            "NRS2",
            None,
            {
                "lrange": 50,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "require_ngood": 15,
                "spline_bkpt": 68,
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "NRS_SLIT",
            "NRS2",
            None,
            {
                "lrange": 50,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "require_ngood": 15,
                "spline_bkpt": 68,
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "NRS_MOS",
            "NRS1",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "require_ngood": 8,
                "spline_bkpt": 30,
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "MIR_MRS",
            "MIRIFUSHORT",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 9, 8, 7, 6],
                "require_ngood": 8,
                "spline_bkpt": 36,
                "space_ratio": 1.2,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "MIR_MRS",
            "MIRIFULONG",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 9, 8, 7, 6],
                "require_ngood": 8,
                "spline_bkpt": 36,
                "space_ratio": 1.2,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
            },
        ),
        (
            "MIR_LRS_SLIT",
            "MIRIMAGE",
            None,
            {
                "lrange": 5,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "require_ngood": 8,
                "spline_bkpt": 40,
                "space_ratio": 1.2,
                "sigma_low": 3.0,
                "sigma_high": 3.0,
                "fit_iter": 3,
            },
        ),
        (
            "MIR_LRS_SLITLESS",
            "MIRIMAGE",
            None,
            {
                "lrange": 5,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "require_ngood": 8,
                "spline_bkpt": 60,
                "space_ratio": 1.2,
                "sigma_low": 3.0,
                "sigma_high": 3.0,
                "fit_iter": 2,
            },
        ),
    ],
)
def test_set_fit_kwargs(mode, detector, slit, expected):
    fit_kwargs = tm._set_fit_kwargs(mode, detector, slit, 10)
    assert list(fit_kwargs.keys()) == list(expected.keys())
    for key in expected:
        if key == "col_index":
            np.testing.assert_equal(list(fit_kwargs[key]), expected[key])
        else:
            assert fit_kwargs[key] == expected[key]


def test_set_fit_kwargs_error():
    with pytest.raises(ValueError, match="Unknown detector"):
        tm._set_fit_kwargs("NRS_SLIT", "NIS", None, 10)


@pytest.mark.parametrize(
    "mode, detector,expected",
    [
        (
            "NRS_IFU",
            "NRS1",
            {
                "pad": 2,
                "trim_ends": True,
            },
        ),
        (
            "NRS_SLIT",
            "NRS1",
            {
                "pad": 1,
                "trim_ends": True,
            },
        ),
        (
            "MIR_MRS",
            "MIRIFUSHORT",
            {
                "pad": 3,
                "trim_ends": False,
            },
        ),
    ],
)
def test_set_oversample_kwargs(mode, detector, expected):
    oversample_kwargs = tm._set_oversample_kwargs(mode, detector)
    assert oversample_kwargs == expected


def test_set_oversample_kwargs_error():
    with pytest.raises(ValueError, match="Unknown detector"):
        tm._set_oversample_kwargs("NRS_SLIT", "NIS")


def test_fit_and_oversample_unsupported_mode():
    model = ImageModel((10, 10))
    model.meta.instrument.detector = "NIS"
    with pytest.raises(ValueError, match="Unknown detector"):
        tm.fit_and_oversample(model)
