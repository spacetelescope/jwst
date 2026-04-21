import numpy as np
from numpy.testing import assert_allclose

from jwst.wfss_contam.disperse import disperse

_SENS_WAVES = np.linspace(1.708, 2.28, 100)
_WMIN, _WMAX = _SENS_WAVES[0], _SENS_WAVES[-1]
_NAXIS = (300, 500)
_SOURCE_ID = 50


def _disperse_one_pixel(grism_wcs, direct_image_wcs, flux_models=None):
    """Run disperse() for a single pixel with flat unit sensitivity."""
    xs = np.array([200])
    ys = np.array([200])
    fluxes = np.array([1.0])
    source_ids = np.array([_SOURCE_ID])
    sens_resp = np.ones_like(_SENS_WAVES)
    return disperse(
        xs,
        ys,
        fluxes,
        source_ids,
        1,
        _WMIN,
        _WMAX,
        _SENS_WAVES,
        sens_resp,
        direct_image_wcs,
        grism_wcs,
        _NAXIS,
        flux_models=flux_models,
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
    assert "model_counts" not in src[source_id[0]]


def test_flux_models_output_structure(grism_wcs, direct_image_with_gradient):
    """Each source dict should have a 'model_counts' list with one entry per model."""
    models = [lambda x: x, lambda x: x**2]
    src = _disperse_one_pixel(grism_wcs, direct_image_with_gradient.meta.wcs, flux_models=models)
    assert "model_counts" in src[_SOURCE_ID]
    assert len(src[_SOURCE_ID]["model_counts"]) == 2
    for mc in src[_SOURCE_ID]["model_counts"]:
        assert mc.shape == src[_SOURCE_ID]["image"].shape


def test_flux_models_scalar_scaling(grism_wcs, direct_image_with_gradient):
    """Scaling a flux model by a constant should scale the model_counts image by the same factor."""
    wcs = direct_image_with_gradient.meta.wcs
    src1 = _disperse_one_pixel(grism_wcs, wcs, flux_models=[lambda x: x])
    src3 = _disperse_one_pixel(grism_wcs, wcs, flux_models=[lambda x: 3.0 * x])
    assert_allclose(src3[_SOURCE_ID]["model_counts"][0], 3.0 * src1[_SOURCE_ID]["model_counts"][0])


def test_flux_models_superposition(grism_wcs, direct_image_with_gradient):
    """The sum of two separate model images should equal the image from their combined model."""
    wcs = direct_image_with_gradient.meta.wcs
    f1 = lambda x: x
    f2 = lambda x: x**2

    src_sep = _disperse_one_pixel(grism_wcs, wcs, flux_models=[f1, f2])
    src_sum = _disperse_one_pixel(grism_wcs, wcs, flux_models=[lambda x: f1(x) + f2(x)])

    combined = src_sep[_SOURCE_ID]["model_counts"][0] + src_sep[_SOURCE_ID]["model_counts"][1]
    assert_allclose(combined, src_sum[_SOURCE_ID]["model_counts"][0], rtol=1e-10)
