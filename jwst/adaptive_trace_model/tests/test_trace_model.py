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
def fit_region_input_no_source_noisy_edge():
    with helpers.nirspec_slit_model() as model:
        flux = model.slits[0].data
        wcs = model.slits[0].meta.wcs
        region_map = (~np.isnan(model.slits[0].wavelength)).astype(int)

    x, y = np.meshgrid(np.arange(flux.shape[1]), np.arange(flux.shape[0]))
    _, alpha, _ = wcs.transform("detector", "slit_frame", x, y)
    err = flux / 2  # for SNR=2 everywhere

    # Add noisy edges
    flux[alpha < -0.5] = 15
    flux[alpha > 0.5] = -15
    err[alpha < -0.5] = 30
    err[alpha > 0.5] = 30

    return flux, err, alpha, region_map


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


def test_crossdisp_profile(fit_2d_spline_input):
    slit, flux, alpha = fit_2d_spline_input
    err = flux * 0.01

    alpha_xdisp, flux_xdisp, err_xdisp, snr_xdisp = tm._crossdisp_profile(flux, err, alpha)
    np.testing.assert_allclose(alpha_xdisp.min(), np.nanmin(alpha))
    np.testing.assert_allclose(alpha_xdisp.max(), np.nanmax(alpha), atol=0.016)
    for xdisp in [flux_xdisp, err_xdisp, snr_xdisp]:
        assert xdisp.shape == alpha_xdisp.shape
    np.testing.assert_allclose(snr_xdisp, 100)


def test_crossdisp_profile_no_error(fit_2d_spline_input):
    slit, flux, alpha = fit_2d_spline_input
    err = np.zeros_like(flux)

    profile = tm._crossdisp_profile(flux, err, alpha)

    # Profile is all None if no errors provided
    for xdisp in profile:
        assert xdisp is None


def test_trim_edges(fit_2d_spline_input):
    slit, flux, alpha = fit_2d_spline_input
    err = flux * 0.01
    alpha_xdisp, _, err_xdisp, snr_xdisp = tm._crossdisp_profile(flux, err, alpha)

    flux_copy = flux.copy()
    tm._trim_edges(flux_copy, alpha, alpha_xdisp, err_xdisp, snr_xdisp)

    # no change for clean edges
    np.testing.assert_allclose(flux_copy, flux)

    # mark edges bad: SNR low, error high
    err_xdisp[:4] = 1000
    err_xdisp[-4:] = 1000
    snr_xdisp[:4] = 4
    snr_xdisp[-4:] = 4
    tm._trim_edges(flux_copy, alpha, alpha_xdisp, err_xdisp, snr_xdisp)

    # each column now has a few NaNs at top and bottom edge of slit
    assert np.sum(np.isnan(flux)) == 0
    assert np.allclose(np.sum(np.isnan(flux_copy), axis=0), 3, atol=1)


def test_trim_edges_all_invalid(fit_2d_spline_input):
    slit, flux, alpha = fit_2d_spline_input
    err = flux * 0.01
    alpha_xdisp, _, err_xdisp, snr_xdisp = tm._crossdisp_profile(flux, err, alpha)

    # all invalid snr
    snr_xdisp[:] = np.nan
    flux_copy = flux.copy()
    tm._trim_edges(flux_copy, alpha, alpha_xdisp, err_xdisp, snr_xdisp)

    # no change
    np.testing.assert_allclose(flux_copy, flux)


def test_threshold_test():
    flux_xdisp = [1, 2, 3, np.nan]
    snr_xdisp = [4, 5, 6, np.nan]
    regnum = 1
    peak_threshold = {1: 2}
    snr_thresold = 5
    assert tm._threshold_test(flux_xdisp, snr_xdisp, regnum, None, None)
    assert tm._threshold_test(flux_xdisp, snr_xdisp, regnum, peak_threshold, None)
    assert tm._threshold_test(flux_xdisp, snr_xdisp, regnum, None, snr_thresold)
    assert tm._threshold_test(flux_xdisp, snr_xdisp, regnum, peak_threshold, snr_thresold)
    assert not tm._threshold_test(flux_xdisp, snr_xdisp, regnum, None, 6)
    assert not tm._threshold_test(flux_xdisp, snr_xdisp, regnum, {1: 3}, None)
    assert not tm._threshold_test(flux_xdisp, snr_xdisp, regnum, {1: 3}, 6)


def test_fit_one_region_below_threshold(caplog, fit_region_input_no_source_noisy_edge):
    flux, err, alpha, region_map = fit_region_input_no_source_noisy_edge
    regnum = 1
    peak_threshold = None

    # SNR threshold below all signal: expect no splines fit
    snr_threshold = 3
    splines = tm._fit_one_region(
        flux, err, alpha, region_map, regnum, peak_threshold, snr_threshold
    )
    assert len(splines) == 0

    # Peak threshold above edge noise but below real signal: expect no splines fit,
    # since the edges are trimmed
    peak_threshold = {1: 10}
    snr_threshold = None
    splines = tm._fit_one_region(
        flux, err, alpha, region_map, regnum, peak_threshold, snr_threshold
    )
    assert len(splines) == 0


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
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": None,
                "require_ngood": None,
                "auto_bkpt_factor": 2.0,
                "auto_ngood_factor": 0.5,
            },
        ),
        (
            "NRS_IFU",
            "NRS2",
            None,
            {
                "lrange": 50,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": None,
                "require_ngood": None,
                "auto_bkpt_factor": 2.0,
                "auto_ngood_factor": 0.5,
            },
        ),
        (
            "NRS_SLIT",
            "NRS2",
            None,
            {
                "lrange": 50,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": None,
                "require_ngood": None,
                "auto_bkpt_factor": 2.0,
                "auto_ngood_factor": 0.5,
            },
        ),
        (
            "NRS_SLIT",
            "NRS1",
            "PRISM",
            {
                "lrange": 10,
                "col_index": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": None,
                "require_ngood": None,
                "auto_bkpt_factor": 1.0,
                "auto_ngood_factor": 0.25,
            },
        ),
        (
            "NRS_MOS",
            "NRS1",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "space_ratio": 1.6,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": None,
                "require_ngood": None,
                "auto_bkpt_factor": 2.0,
                "auto_ngood_factor": 0.5,
            },
        ),
        (
            "MIR_MRS",
            "MIRIFUSHORT",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 9, 8, 7, 6],
                "space_ratio": 1.2,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": 36,
                "require_ngood": 8,
                "auto_bkpt_factor": None,
                "auto_ngood_factor": None,
            },
        ),
        (
            "MIR_MRS",
            "MIRIFULONG",
            None,
            {
                "lrange": 50,
                "col_index": [0, 1, 2, 3, 4, 5, 9, 8, 7, 6],
                "space_ratio": 1.2,
                "sigma_low": 2.5,
                "sigma_high": 2.5,
                "fit_iter": 3,
                "spline_bkpt": 36,
                "require_ngood": 8,
                "auto_bkpt_factor": None,
                "auto_ngood_factor": None,
            },
        ),
        (
            "MIR_LRS_SLIT",
            "MIRIMAGE",
            None,
            {
                "lrange": 5,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "space_ratio": 1.2,
                "sigma_low": 3.0,
                "sigma_high": 3.0,
                "fit_iter": 3,
                "spline_bkpt": 40,
                "require_ngood": 8,
                "auto_bkpt_factor": None,
                "auto_ngood_factor": None,
            },
        ),
        (
            "MIR_LRS_SLITLESS",
            "MIRIMAGE",
            None,
            {
                "lrange": 5,
                "col_index": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                "space_ratio": 1.2,
                "sigma_low": 3.0,
                "sigma_high": 3.0,
                "fit_iter": 2,
                "spline_bkpt": 60,
                "require_ngood": 8,
                "auto_bkpt_factor": None,
                "auto_ngood_factor": None,
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
            {"pad": 2, "trim_ends": True, "require_ngood": 8},
        ),
        (
            "NRS_SLIT",
            "NRS1",
            {"pad": 1, "trim_ends": True, "require_ngood": 8},
        ),
        (
            "MIR_MRS",
            "MIRIFUSHORT",
            {"pad": 3, "trim_ends": False, "require_ngood": 8},
        ),
    ],
)
def test_set_oversample_kwargs(mode, detector, expected):
    oversample_kwargs = tm._set_oversample_kwargs(mode, detector)
    assert oversample_kwargs == expected


def test_set_oversample_kwargs_error():
    with pytest.raises(ValueError, match="Unknown detector"):
        tm._set_oversample_kwargs("NRS_SLIT", "NIS")


def test_fit_and_oversample_unknown_detector():
    model = ImageModel((10, 10))
    model.meta.instrument.detector = "NIS"
    with pytest.raises(ValueError, match="Unknown detector"):
        tm.fit_and_oversample(model)


def test_fit_and_oversample_unsupported_exptype():
    model = ImageModel((10, 10))
    model.meta.exposure.type = "MIR_WFSS"
    model.meta.instrument.detector = "MIRIMAGE"
    with pytest.raises(ValueError, match="MIRI WFSS is not supported"):
        tm.fit_and_oversample(model)
