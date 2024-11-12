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
    profile[3] = 0.5    # lower limit: middle of pixel 3
    profile[4:7] = 1.0  # center of aperture
    profile[7] = 0.5    # upper limit: middle of pixel 7

    return (image, var_rnoise, var_poisson, var_rflat,
            profile, weights, profile_bg)


def test_extract_src_flux(inputs_constant):
    (image, var_rnoise, var_poisson, var_rflat,
     profile, weights, profile_bg) = inputs_constant

    (total_flux, f_var_rnoise, f_var_poisson, f_var_flat,
     bkg_flux, b_var_rnoise, b_var_poisson, b_var_flat,
     npixels, model) = extract1d.extract1d(
        image, [profile], var_rnoise, var_poisson, var_rflat,
        weights=weights, profile_bg=profile_bg)

    # check the value at column 2
    # 0.5 * 17. + 22. + 27. + 32. + 0.5 * 37.
    assert math.isclose(total_flux[0][2], 108., rel_tol=1.e-8, abs_tol=1.e-8)
    assert bkg_flux[0][2] == 0.
    assert math.isclose(npixels[0][2], 4., rel_tol=1.e-8, abs_tol=1.e-8)

    # set a NaN value in the column of interest
    image[5, 2] = np.nan

    (total_flux, f_var_rnoise, f_var_poisson, f_var_flat,
     bkg_flux, b_var_rnoise, b_var_poisson, b_var_flat,
     npixels, model) = extract1d.extract1d(
        image, [profile], var_rnoise, var_poisson, var_rflat,
        weights=weights, profile_bg=profile_bg)

    # 0.5 * 17. + 22. + 32. + 0.5 * 37.
    assert math.isclose(total_flux[0][2], 81., rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(npixels[0][2], 3., rel_tol=1.e-8, abs_tol=1.e-8)

    # set the whole column to NaN
    image[:, 2] = np.nan

    (total_flux, f_var_rnoise, f_var_poisson, f_var_flat,
     bkg_flux, b_var_rnoise, b_var_poisson, b_var_flat,
     npixels, model) = extract1d.extract1d(
        image, [profile], var_rnoise, var_poisson, var_rflat,
        weights=weights, profile_bg=profile_bg)

    assert np.isnan(total_flux[0][2])
    assert npixels[0][2] == 0.


def test_extract_src_flux_empty_interval(inputs_constant):
    (image, var_rnoise, var_poisson, var_rflat,
     profile, weights, profile_bg) = inputs_constant

    # empty extraction range
    profile[:] = 0.0

    (total_flux, f_var_rnoise, f_var_poisson, f_var_flat,
     bkg_flux, b_var_rnoise, b_var_poisson, b_var_flat,
     npixels, model) = extract1d.extract1d(
        image, [profile], var_rnoise, var_poisson, var_rflat,
        weights=weights, profile_bg=profile_bg)

    # empty interval, so no flux returned
    assert np.all(np.isnan(total_flux))
    assert np.all(bkg_flux == 0.)
    assert np.all(npixels == 0.)
