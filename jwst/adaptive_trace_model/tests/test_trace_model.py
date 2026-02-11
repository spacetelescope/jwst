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
    splines, scales = tm.fit_2d_spline_trace(flux, alpha)

    assert len(splines) == flux.shape[1]
    assert len(scales) == flux.shape[1]

    # for unscaled fits, output scales should be close to 1
    np.testing.assert_allclose(list(scales.values()), 1.0, atol=0.06)

    # evaluated fits should be close to input data
    region_map = (~np.isnan(slit.wavelength)).astype(int)
    trace_used, full_trace = tm._trace_image(
        flux.shape, {1: splines}, {1: scales}, region_map, alpha
    )

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
    splines, scales = tm.fit_2d_spline_trace(flux, alpha)
    assert len(splines) == 0
    assert len(scales) == 0
    assert "Spline fit failed" in caplog.text


@pytest.mark.parametrize(
    "detector,expected",
    [
        (
            "NRS1",
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "require_ngood": 15,
                "spline_bkpt": 62,
                "space_ratio": 1.6,
            },
        ),
        (
            "NRS2",
            {
                "lrange": 50,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "require_ngood": 15,
                "spline_bkpt": 62,
                "space_ratio": 1.6,
            },
        ),
        (
            "MIRIFUSHORT",
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 9, 8, 7, 6],
                "require_ngood": 8,
                "spline_bkpt": 36,
                "space_ratio": 1.2,
            },
        ),
        (
            "MIRIFULONG",
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 9, 8, 7, 6],
                "require_ngood": 8,
                "spline_bkpt": 36,
                "space_ratio": 1.2,
            },
        ),
    ],
)
def test_set_fit_kwargs(detector, expected):
    fit_kwargs = tm._set_fit_kwargs(detector, 10)
    assert list(fit_kwargs.keys()) == list(expected.keys())
    for key in expected:
        if key == "col_index":
            np.testing.assert_equal(list(fit_kwargs[key]), expected[key])
        else:
            assert fit_kwargs[key] == expected[key]


def test_set_fit_kwargs_error():
    with pytest.raises(ValueError, match="Unknown detector"):
        tm._set_fit_kwargs("NIS", 10)


@pytest.mark.parametrize(
    "detector,expected",
    [
        (
            "NRS1",
            {
                "pad": 2,
                "trim_ends": True,
            },
        ),
        (
            "MIRIFUSHORT",
            {
                "pad": 3,
                "trim_ends": False,
            },
        ),
    ],
)
def test_set_oversample_kwargs(detector, expected):
    oversample_kwargs = tm._set_oversample_kwargs(detector)
    assert oversample_kwargs == expected


def test_set_oversample_kwargs_error():
    with pytest.raises(ValueError, match="Unknown detector"):
        tm._set_oversample_kwargs("NIS")


@pytest.mark.parametrize(
    "detector,message",
    [("NRS1", "Unsupported mode"), ("MIRIFULONG", "Unsupported mode"), ("NIS", "Unknown detector")],
)
def test_fit_and_oversample_unsupported_model(detector, message):
    model = ImageModel()
    model.meta.instrument.detector = detector
    with pytest.raises(ValueError, match=message):
        tm.fit_and_oversample(model)
