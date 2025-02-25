"""
Test for extract_1d._extract_src_flux
"""

import math

import numpy as np
import pytest

from jwst.extract_1d import extract1d


@pytest.fixture
def inputs_constant():
    shape = (9, 5)
    image = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    var_rnoise = image.copy()
    var_poisson = image.copy()
    var_rflat = image.copy()
    weights = None
    profile_bg = None

    profile = np.zeros_like(image)
    profile[3] = 0.5  # lower limit: middle of pixel 3
    profile[4:7] = 1.0  # center of aperture
    profile[7] = 0.5  # upper limit: middle of pixel 7

    return (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg)


@pytest.fixture
def inputs_with_source():
    shape = (9, 5)
    image = np.full(shape, 0.0)
    image[3] = 5.0
    image[4] = 10.0
    image[5] = 5.0
    var_rnoise = np.full(shape, 0.1)
    var_poisson = image * 0.05
    var_rflat = image * 0.05
    weights = None
    profile_bg = None

    profile = np.zeros_like(image)
    profile[3] = 0.25
    profile[4] = 0.5
    profile[5] = 0.25

    return (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg)


def test_extract_src_flux(inputs_constant):
    (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg) = inputs_constant

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
        image, [profile], var_rnoise, var_poisson, var_rflat, weights=weights, profile_bg=profile_bg
    )

    # check the value at column 2
    # 0.5 * 17. + 22. + 27. + 32. + 0.5 * 37.
    assert math.isclose(total_flux[0][2], 108.0, rel_tol=1.0e-8, abs_tol=1.0e-8)
    assert bkg_flux[0][2] == 0.0
    assert math.isclose(npixels[0][2], 4.0, rel_tol=1.0e-8, abs_tol=1.0e-8)

    # set a NaN value in the column of interest
    image[5, 2] = np.nan

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
        image, [profile], var_rnoise, var_poisson, var_rflat, weights=weights, profile_bg=profile_bg
    )

    # 0.5 * 17. + 22. + 32. + 0.5 * 37.
    assert math.isclose(total_flux[0][2], 81.0, rel_tol=1.0e-8, abs_tol=1.0e-8)
    assert math.isclose(npixels[0][2], 3.0, rel_tol=1.0e-8, abs_tol=1.0e-8)

    # set the whole column to NaN
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
        image, [profile], var_rnoise, var_poisson, var_rflat, weights=weights, profile_bg=profile_bg
    )

    assert np.isnan(total_flux[0][2])
    assert npixels[0][2] == 0.0


def test_extract_src_flux_empty_interval(inputs_constant):
    (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg) = inputs_constant

    # empty extraction range
    profile[:] = 0.0

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
        image, [profile], var_rnoise, var_poisson, var_rflat, weights=weights, profile_bg=profile_bg
    )

    # empty interval, so no flux returned
    assert np.all(np.isnan(total_flux))
    assert np.all(bkg_flux == 0.0)
    assert np.all(npixels == 0.0)


@pytest.mark.parametrize("use_weights", [True, False])
def test_extract_optimal(inputs_with_source, use_weights):
    (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg) = inputs_with_source

    if use_weights:
        weights = 1 / var_rnoise

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
        extraction_type="optimal",
    )

    # Total flux should be well modeled
    assert np.allclose(total_flux[0], 20)
    assert np.allclose(bkg_flux[0], 0)
    assert np.allclose(npixels[0], 2.66667)

    # set a NaN value in a column of interest
    image[4, 2] = np.nan
    if use_weights:
        weights[4, 2] = 0

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
        extraction_type="optimal",
    )

    # Total flux is still well modeled from 2 pixels
    assert np.allclose(total_flux[0], 20)
    assert np.isclose(npixels[0, 2], 1.33333)

    # set the whole column to NaN
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
        image, [profile], var_rnoise, var_poisson, var_rflat, weights=weights, profile_bg=profile_bg
    )

    # Now the flux can no longer be estimated in that column
    assert np.isnan(total_flux[0][2])
    assert npixels[0][2] == 0.0


def test_too_many_profiles(inputs_constant):
    (image, var_rnoise, var_poisson, var_rflat, profile, weights, profile_bg) = inputs_constant

    with pytest.raises(ValueError, match="not supported with 2 input profiles"):
        extract1d.extract1d(image, [profile, profile.copy()], var_rnoise, var_poisson, var_rflat)
