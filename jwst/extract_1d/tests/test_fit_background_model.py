"""
Tests for extract_1d background fitting
"""

import numpy as np
import pytest

from jwst.extract_1d import extract1d


@pytest.fixture
def inputs_constant():
    shape = (9, 5)
    image = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    var_poisson = image.copy()
    var_rnoise = image.copy()
    var_rflat = image.copy()

    profile = np.zeros_like(image)
    profile[1] = 1.0  # one pixel aperture
    weights = None

    profile_bg = np.zeros_like(image)
    profile_bg[4:8] = 1.0  # background region
    bkg_order = 0
    bkg_fit = "poly"

    return (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    )


@pytest.fixture
def inputs_with_source():
    shape = (9, 5)
    image = np.full(shape, 1.0)
    image[3] += 5.0
    image[4] += 10.0
    image[5] += 5.0

    var_rnoise = np.full(shape, 0.1)
    var_poisson = image * 0.05
    var_rflat = image * 0.05
    weights = 1 / var_rnoise

    # Most of the image is set to a low but non-zero weight
    # (contribution to PSF is small)
    profile = np.full_like(image, 1e-6)

    # The source crosses 3 pixels, with the central pixel
    # twice as high as the other two
    profile[3] = 1
    profile[4] = 2
    profile[5] = 1

    # Normalize across the spatial dimension
    profile /= np.sum(profile, axis=0)

    profile_bg = None

    return (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg)


@pytest.mark.parametrize("bkg_fit_type", ["poly", "median"])
def test_fit_background(inputs_constant, bkg_fit_type):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant

    (
        total_flux,
        f_var_rnoise,
        f_var_poisson,
        f_var_flat,
        bkg_flux,
        b_var_rnoise,
        b_var_poisson,
        b_var_flat,
        npixels,
        model,
    ) = extract1d.extract1d(
        image,
        [profile],
        var_rnoise,
        var_poisson,
        var_rflat,
        weights=weights,
        profile_bg=profile_bg,
        fit_bkg=True,
        bkg_fit_type=bkg_fit_type,
        bkg_order=bkg_order,
    )

    if bkg_fit_type == "median":
        extra_factor = 1.2**2
    else:
        extra_factor = 1.0

    assert np.allclose(bkg_flux[0], np.sum(image[4:8], axis=0) / 4)
    assert np.allclose(b_var_rnoise[0], extra_factor * np.sum(image[4:8], axis=0) / 16)
    assert np.allclose(b_var_poisson[0], extra_factor * np.sum(image[4:8], axis=0) / 16)
    assert np.allclose(b_var_flat[0], extra_factor * np.sum(image[4:8], axis=0) / 16)


@pytest.mark.parametrize("bkg_fit_type", ["poly", "median"])
def test_fit_background_with_smoothing(inputs_constant, bkg_fit_type):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant

    (
        total_flux,
        f_var_rnoise,
        f_var_poisson,
        f_var_flat,
        bkg_flux,
        b_var_rnoise,
        b_var_poisson,
        b_var_flat,
        npixels,
        model,
    ) = extract1d.extract1d(
        image,
        [profile],
        var_rnoise,
        var_poisson,
        var_rflat,
        weights=weights,
        profile_bg=profile_bg,
        fit_bkg=True,
        bkg_fit_type=bkg_fit_type,
        bkg_order=bkg_order,
        bg_smooth_length=3,
    )

    if bkg_fit_type == "median":
        extra_factor = 1.2**2
    else:
        extra_factor = 1.0

    # Should be the same as background fit without smoothing,
    # except for edge effects
    expected = np.sum(image[4:8], axis=0) / 4
    assert np.allclose(bkg_flux[0][1:-1], expected[1:-1])
    assert np.allclose(b_var_rnoise[0][1:-1], extra_factor * expected[1:-1] / 4)
    assert np.allclose(b_var_poisson[0][1:-1], extra_factor * expected[1:-1] / 4)
    assert np.allclose(b_var_flat[0][1:-1], extra_factor * expected[1:-1] / 4)


def test_handles_nan(inputs_constant):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant
    image[:, 2] = np.nan

    (
        total_flux,
        f_var_rnoise,
        f_var_poisson,
        f_var_flat,
        bkg_flux,
        b_var_rnoise,
        b_var_poisson,
        b_var_flat,
        npixels,
        model,
    ) = extract1d.extract1d(
        image,
        [profile],
        var_rnoise,
        var_poisson,
        var_rflat,
        weights=weights,
        profile_bg=profile_bg,
        fit_bkg=True,
        bkg_fit_type=bkg_fit,
        bkg_order=bkg_order,
    )

    assert np.allclose(bkg_flux[0], np.nansum(image[4:8], axis=0) / 4)
    assert np.allclose(b_var_rnoise[0], np.nansum(image[4:8], axis=0) / 16)
    assert np.allclose(b_var_poisson[0], np.nansum(image[4:8], axis=0) / 16)
    assert np.allclose(b_var_flat[0], np.nansum(image[4:8], axis=0) / 16)


@pytest.mark.parametrize("bkg_order_val", [0, 1, 2])
def test_handles_one_value(inputs_constant, bkg_order_val):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant
    profile_bg[:] = 0.0
    profile_bg[4:6] = 1.0
    image[4] = np.nan

    (
        total_flux,
        f_var_rnoise,
        f_var_poisson,
        f_var_flat,
        bkg_flux,
        b_var_rnoise,
        b_var_poisson,
        b_var_flat,
        npixels,
        model,
    ) = extract1d.extract1d(
        image,
        [profile],
        var_rnoise,
        var_poisson,
        var_rflat,
        weights=weights,
        profile_bg=profile_bg,
        fit_bkg=True,
        bkg_fit_type=bkg_fit,
        bkg_order=bkg_order_val,
    )

    if bkg_order_val == 0:
        expected = image[5]
    else:
        # not enough input for fit: value is set to 0.0
        expected = 0.0

    assert np.allclose(bkg_flux[0], expected)
    assert np.allclose(b_var_rnoise[0], expected)
    assert np.allclose(b_var_poisson[0], expected)
    assert np.allclose(b_var_flat[0], expected)


def test_handles_empty_interval(inputs_constant):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant
    profile_bg[:] = 0.0

    (
        total_flux,
        f_var_rnoise,
        f_var_poisson,
        f_var_flat,
        bkg_flux,
        b_var_rnoise,
        b_var_poisson,
        b_var_flat,
        npixels,
        model,
    ) = extract1d.extract1d(
        image,
        [profile],
        var_rnoise,
        var_poisson,
        var_rflat,
        weights=weights,
        profile_bg=profile_bg,
        fit_bkg=True,
        bkg_fit_type=bkg_fit,
        bkg_order=bkg_order,
    )

    assert np.allclose(bkg_flux[0], 0.0)
    assert np.allclose(b_var_rnoise[0], 0.0)
    assert np.allclose(b_var_poisson[0], 0.0)
    assert np.allclose(b_var_flat[0], 0.0)


def test_bad_fit_type(inputs_constant):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant
    with pytest.raises(ValueError, match="bkg_fit_type should be 'median' or 'poly'"):
        extract1d.extract1d(
            image,
            [profile],
            var_rnoise,
            var_poisson,
            var_rflat,
            weights=weights,
            profile_bg=profile_bg,
            fit_bkg=True,
            bkg_fit_type="mean",
            bkg_order=bkg_order,
        )


@pytest.mark.parametrize("smooth", [1.5, 2])
def test_smooth_length(inputs_constant, smooth):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant
    with pytest.raises(ValueError, match="should be an odd integer >= 1"):
        extract1d.extract1d(
            image,
            [profile],
            var_rnoise,
            var_poisson,
            var_rflat,
            weights=weights,
            profile_bg=profile_bg,
            fit_bkg=True,
            bkg_fit_type=bkg_fit,
            bkg_order=bkg_order,
            bg_smooth_length=smooth,
        )


@pytest.mark.parametrize("extraction_type", ["box", "optimal"])
@pytest.mark.parametrize("bkg_order_val", [-1, 2.3])
def test_bad_fit_order(inputs_constant, extraction_type, bkg_order_val):
    (
        image,
        var_rnoise,
        var_poisson,
        var_rflat,
        profile,
        weights,
        profile_bg,
        bkg_fit,
        bkg_order,
    ) = inputs_constant
    with pytest.raises(ValueError, match="bkg_order must be an integer >= 0"):
        extract1d.extract1d(
            image,
            [profile],
            var_rnoise,
            var_poisson,
            var_rflat,
            weights=weights,
            profile_bg=profile_bg,
            fit_bkg=True,
            bkg_fit_type="poly",
            bkg_order=bkg_order_val,
            extraction_type=extraction_type,
        )


@pytest.mark.parametrize("use_weights", [True, False])
@pytest.mark.parametrize("bkg_order_val", [0, 1, 2])
def test_fit_background_optimal(inputs_with_source, use_weights, bkg_order_val):
    (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg) = inputs_with_source

    if not use_weights:
        weights = None

    result = extract1d.extract1d(
        image,
        [profile],
        var_rnoise,
        var_poisson,
        var_rflat,
        weights=weights,
        profile_bg=profile_bg,
        fit_bkg=True,
        bkg_fit_type="poly",
        bkg_order=bkg_order_val,
        extraction_type="optimal",
    )

    flux = result[0][0]
    background = result[4][0]

    assert np.allclose(flux, 20.0)
    assert np.allclose(background, 2.66667)
