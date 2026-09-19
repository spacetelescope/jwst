"""Tests for jwst.wfss_contam.trace_pdt."""

import numpy as np
import pytest
from astropy.modeling.mappings import Mapping
from numpy.testing import assert_allclose

from jwst.wfss_contam.trace_pdt import TracePDT, build_trace_pdt

_ORDER = 1
_WMIN, _WMAX = 1.708, 2.28
_NAXIS = (2048, 2048)


def _exact_transform(grism_wcs):
    """Build the exact "detector" -> "grism_detector" (x, y) transform for comparison."""
    transform = grism_wcs.get_transform("detector", "grism_detector")
    n_outputs = len(transform.outputs)
    return transform | Mapping((0, 1), n_inputs=n_outputs)


def test_trace_pdt_matches_exact_transform(grism_wcs):
    """Interpolated trace positions should closely match the exact transform off-grid."""
    trace_pdt = build_trace_pdt(grism_wcs, _ORDER, _WMIN, _WMAX, _NAXIS, spacing=100)
    assert isinstance(trace_pdt, TracePDT)

    # Sample random points strictly inside the grid bounds, away from grid nodes.
    rng = np.random.default_rng(42)
    n_pts = 200
    x0 = rng.uniform(10, _NAXIS[0] - 10, n_pts)
    y0 = rng.uniform(10, _NAXIS[1] - 10, n_pts)
    wavelength = rng.uniform(_WMIN, _WMAX, n_pts)

    exact_transform = _exact_transform(grism_wcs)
    x_exact, y_exact = exact_transform(x0, y0, wavelength, _ORDER)
    x_interp, y_interp = trace_pdt(x0, y0, wavelength)

    assert_allclose(x_interp, x_exact, atol=0.05)
    assert_allclose(y_interp, y_exact, atol=0.05)


def test_trace_pdt_preserves_shape(grism_wcs):
    """The TracePDT should return outputs with the same shape as the inputs."""
    trace_pdt = build_trace_pdt(grism_wcs, _ORDER, _WMIN, _WMAX, _NAXIS, spacing=100)

    x0 = np.full((5, 7), 500.0)
    y0 = np.full((5, 7), 600.0)
    wavelength = np.full((5, 7), 2.0)
    x, y = trace_pdt(x0, y0, wavelength)
    assert x.shape == (5, 7)
    assert y.shape == (5, 7)


def test_evaluate_grid_matches_call(grism_wcs):
    """evaluate_grid should agree with the generic __call__ on the same query points."""
    trace_pdt = build_trace_pdt(grism_wcs, _ORDER, _WMIN, _WMAX, _NAXIS, spacing=100)

    rng = np.random.default_rng(0)
    n_pixels = 37
    n_wave = 13
    x0 = rng.uniform(10, _NAXIS[0] - 10, n_pixels)
    y0 = rng.uniform(10, _NAXIS[1] - 10, n_pixels)
    wavelength = rng.uniform(_WMIN, _WMAX, n_wave)

    x_grid, y_grid = trace_pdt.evaluate_grid(x0, y0, wavelength)
    assert x_grid.shape == (n_wave, n_pixels)
    assert y_grid.shape == (n_wave, n_pixels)

    # Build the same (n_wave, n_pixels) outer product manually and compare to __call__.
    x0_rep = np.repeat(x0[np.newaxis, :], n_wave, axis=0)
    y0_rep = np.repeat(y0[np.newaxis, :], n_wave, axis=0)
    lam_rep = np.repeat(wavelength[:, np.newaxis], n_pixels, axis=1)
    x_call, y_call = trace_pdt(x0_rep, y0_rep, lam_rep)

    assert_allclose(x_grid, x_call, rtol=1e-10)
    assert_allclose(y_grid, y_call, rtol=1e-10)


def test_exact_wavelength_grid_matches_coarse_grid(grism_wcs):
    """Passing lam_grid should skip wavelength interpolation but agree with the coarse grid."""
    lam_grid = np.linspace(_WMIN, _WMAX, 50)
    trace_pdt_exact = build_trace_pdt(
        grism_wcs, _ORDER, _WMIN, _WMAX, _NAXIS, spacing=100, lam_grid=lam_grid
    )
    trace_pdt_coarse = build_trace_pdt(grism_wcs, _ORDER, _WMIN, _WMAX, _NAXIS, spacing=100)
    assert trace_pdt_exact._exact_wavelength_grid is True
    assert trace_pdt_coarse._exact_wavelength_grid is False

    x0 = np.array([500.0, 800.0])
    y0 = np.array([600.0, 900.0])

    x_exact, y_exact = trace_pdt_exact.evaluate_grid(x0, y0, lam_grid)
    x_coarse, y_coarse = trace_pdt_coarse.evaluate_grid(x0, y0, lam_grid)

    assert_allclose(x_exact, x_coarse, atol=0.05)
    assert_allclose(y_exact, y_coarse, atol=0.05)


def test_evaluate_grid_exact_wavelength_grid_raises_on_mismatch(grism_wcs):
    """evaluate_grid must not silently interpolate a wavelength array it wasn't built for."""
    lam_grid = np.linspace(_WMIN, _WMAX, 50)
    trace_pdt = build_trace_pdt(
        grism_wcs, _ORDER, _WMIN, _WMAX, _NAXIS, spacing=100, lam_grid=lam_grid
    )

    x0 = np.array([500.0])
    y0 = np.array([600.0])
    wrong_lam = np.linspace(_WMIN, _WMAX, 30)
    with pytest.raises(ValueError, match="exact_wavelength_grid"):
        trace_pdt.evaluate_grid(x0, y0, wrong_lam)
