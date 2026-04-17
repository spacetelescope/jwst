import numpy as np
from numpy.testing import assert_allclose

from jwst.wfss_contam.disperse import disperse


def test_disperse_oversample_same_result(grism_wcs, direct_image_with_gradient):
    """Coverage for bug where wavelength oversampling led to double-counted fluxes."""
    x0 = np.array([200.5])
    y0 = np.array([200.5])
    order = 1
    flxs = np.array([[1.0]])
    band_wave = np.array([2.0])
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
            band_wave,
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


def _disperse_one_pixel(fluxes, band_wavelengths, grism_wcs, direct_image_wcs):
    """Run disperse() for a single pixel with unit sensitivity, returning the output image."""
    x0 = np.array([200.5])
    y0 = np.array([200.5])
    source_id = np.array([50])
    naxis = (300, 500)
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = float(sens_waves[0]), float(sens_waves[-1])
    sens_resp = np.ones(100)
    src = disperse(
        x0,
        y0,
        fluxes,
        band_wavelengths,
        source_id,
        order=1,
        wmin=wmin,
        wmax=wmax,
        sens_waves=sens_waves,
        sens_resp=sens_resp,
        direct_image_wcs=direct_image_wcs,
        grism_wcs=grism_wcs,
        naxis=naxis,
    )
    return src[source_id[0]]["image"]


def test_disperse_n1_matches_n2_constant(grism_wcs, direct_image_with_gradient):
    """
    N=1 and N=2 with the same constant flux value should produce identical output.

    This tests that the single-band fallback path and the interp1d path
    agree when the SED is flat.
    """
    direct_image_wcs = direct_image_with_gradient.meta.wcs
    img_n1 = _disperse_one_pixel(
        fluxes=np.array([[1.0]]),
        band_wavelengths=np.array([2.0]),
        grism_wcs=grism_wcs,
        direct_image_wcs=direct_image_wcs,
    )
    img_n2 = _disperse_one_pixel(
        fluxes=np.array([[1.0], [1.0]]),
        band_wavelengths=np.array([1.708, 2.28]),
        grism_wcs=grism_wcs,
        direct_image_wcs=direct_image_wcs,
    )
    assert_allclose(img_n1, img_n2, rtol=1e-10)


def test_disperse_flux_extrapolation(grism_wcs, direct_image_with_gradient):
    """Test that extrapolation holds the edge flux value outside the band coverage."""
    direct_image_wcs = direct_image_with_gradient.meta.wcs
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = float(sens_waves[0]), float(sens_waves[-1])
    wmid1 = wmin + 0.3 * (wmax - wmin)
    wmid2 = wmax - 0.3 * (wmax - wmin)

    img_extrap = _disperse_one_pixel(
        fluxes=np.array([[0.0], [2.0]]),
        band_wavelengths=np.array([wmid1, wmid2]),
        grism_wcs=grism_wcs,
        direct_image_wcs=direct_image_wcs,
    )
    # manually padded 4-point SED that encodes the same flat extrapolation
    img_explicit = _disperse_one_pixel(
        fluxes=np.array([[0.0], [0.0], [2.0], [2.0]]),
        band_wavelengths=np.array([wmin, wmid1, wmid2, wmax]),
        grism_wcs=grism_wcs,
        direct_image_wcs=direct_image_wcs,
    )
    assert_allclose(img_extrap, img_explicit, rtol=1e-10)


def test_disperse_flux_distribution(grism_wcs, direct_image_with_gradient):
    """
    Test that multi-band input is being nontrivially accounted for.

    Use a sensitivity curve that passes only the lower half of the wavelength range,
    so only the short wavelengths contribute to the output. Then test two SEDs
    with the same average flux but different wavelength dependence.
    """
    direct_image_wcs = direct_image_with_gradient.meta.wcs
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = float(sens_waves[0]), float(sens_waves[-1])

    # Sensitivity passes only the lower half of the wavelength range
    # This allows us to test that the interpolation is actually doing something
    sens_resp = np.zeros(100)
    sens_resp[:50] = 1.0

    def run(fluxes, bw):
        src = disperse(
            np.array([200.5]),
            np.array([200.5]),
            fluxes,
            bw,
            np.array([50]),
            order=1,
            wmin=wmin,
            wmax=wmax,
            sens_waves=sens_waves,
            sens_resp=sens_resp,
            direct_image_wcs=direct_image_wcs,
            grism_wcs=grism_wcs,
            naxis=(300, 500),
        )
        return src[50]["image"].sum()

    bw = np.array([wmin, wmax])
    sum_falling = run(np.array([[2.0], [0.0]]), bw)
    sum_rising = run(np.array([[0.0], [2.0]]), bw)

    assert sum_falling / sum_rising > 2.0
