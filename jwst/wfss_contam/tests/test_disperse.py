import numpy as np
from numpy.testing import assert_allclose

from jwst.wfss_contam.disperse import (
    _build_mean_wavelength_image_of_source,
    _collect_outputs_by_source,
    disperse,
)


def test_disperse_oversample_same_result(grism_wcs, direct_image_with_gradient):
    """Coverage for bug where wavelength oversampling led to double-counted fluxes."""
    x0 = np.array([200.5])
    y0 = np.array([200.5])
    order = 1
    flxs = np.array([1.0])
    source_id = np.array([50])
    naxis = (300, 500)
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)
    direct_image_wcs = direct_image_with_gradient.meta.wcs

    output_images = []
    for os in [2, 3]:
        src = disperse(
            x0,
            y0,
            flxs,
            source_id,
            order,
            wmin,
            wmax,
            sens_waves,
            sens_resp,
            direct_image_wcs,
            grism_wcs,
            naxis,
            oversample_factor=os,
        )
        output_images.append(src[source_id[0]]["image"])

    # different oversampling gives different effects at the ends
    # unsure if this is a bug or not, but the middle should definitely be the same
    assert_allclose(output_images[0][2:-2, :], output_images[1][2:-2], rtol=1e-5)


def test_build_mean_image_of_source_duplicates_unequal_weights():
    """Two inputs at the same location with unequal weights: result is the weighted mean."""
    x = np.array([0, 0])
    y = np.array([0, 0])
    values = np.array([1.0, 3.0])
    areas = np.array([3.0, 1.0])  # first value has 3x the area
    bounds = [0, 0, 0, 0]
    result = _build_mean_wavelength_image_of_source(x, y, values, areas, bounds)
    assert result.shape == (1, 1)
    # expected = (3*1 + 1*3) / (3+1) = 6/4 = 1.5
    assert_allclose(result[0, 0], 1.5)


def test_build_mean_image_of_source_empty_pixels_are_zero():
    """Pixels with no input contribution should be zero, not NaN."""
    x = np.array([0])
    y = np.array([0])
    values = np.array([5.0])
    areas = np.array([1.0])
    bounds = [0, 1, 0, 1]  # 2x2 output, only (0,0) filled
    result = _build_mean_wavelength_image_of_source(x, y, values, areas, bounds)
    assert result.shape == (2, 2)
    assert_allclose(result[0, 0], 5.0)
    assert_allclose(result[0, 1], 0.0)
    assert_allclose(result[1, 0], 0.0)
    assert_allclose(result[1, 1], 0.0)


def test_collect_outputs_by_source_wavelengths_present():
    """_collect_outputs_by_source stores a 'wavelengths' array for each source."""
    rng = np.random.default_rng(42)
    n = 20
    xs = rng.integers(0, 5, n)
    ys = rng.integers(0, 5, n)
    lambdas = rng.uniform(1.7, 2.3, n)
    areas = rng.uniform(0.1, 1.0, n)
    counts = rng.uniform(0.0, 1.0, n)
    source_ids = np.ones(n, dtype=int)  # single source

    out = _collect_outputs_by_source(xs, ys, lambdas, areas, counts, source_ids)
    assert 1 in out
    assert "wavelengths" in out[1]
    lam = out[1]["wavelengths"]
    img = out[1]["image"]
    assert lam.shape == img.shape
    # wavelengths should be within the input lambda range or zero (empty pixels)
    assert np.all((lam >= 1.7) | (lam == 0.0))
    assert np.all((lam <= 2.3) | (lam == 0.0))
