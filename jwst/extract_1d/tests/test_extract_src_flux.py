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
    var_poisson = image.copy()
    var_rnoise = image.copy()
    var_rflat = image.copy()
    x = 2
    j = 2
    lam = 1.234                 # an arbitrary value (not actually used)
    weights = None
    bkgmodels = [None, None, None, None]
    lower = np.zeros(shape[1], dtype=np.float64) + 3.   # middle of pixel 3
    upper = np.zeros(shape[1], dtype=np.float64) + 7.   # middle of pixel 7
    srclim = [[lower, upper]]

    return (image, var_poisson, var_rnoise, var_rflat,
            x, j, lam, srclim, weights, bkgmodels)


def test_extract_src_flux(inputs_constant):
    (image, var_poisson, var_rnoise, var_rflat,
     x, j, lam, srclim, weights, bkgmodels) = inputs_constant

    (total_flux, f_var_poisson, f_var_rnoise, f_var_flat,
     bkg_flux, b_var_poisson, b_var_rnoise, b_var_flat,
     tarea, twht) = extract1d._extract_src_flux(
        image, var_poisson, var_rnoise, var_rflat,
        x, j, lam, srclim, weights, bkgmodels)

    # 0.5 * 17. + 22. + 27. + 32. + 0.5 * 37.
    assert math.isclose(total_flux, 108., rel_tol=1.e-8, abs_tol=1.e-8)

    assert bkg_flux == 0.

    assert math.isclose(tarea, 4., rel_tol=1.e-8, abs_tol=1.e-8)

    image[5, 2] = np.nan

    (total_flux, f_var_poisson, f_var_rnoise, f_var_flat,
     bkg_flux, b_var_poisson, b_var_rnoise, b_var_flat,
     tarea, twht) = extract1d._extract_src_flux(
        image, var_poisson, var_rnoise, var_rflat,
        x, j, lam, srclim, weights, bkgmodels)

    # 0.5 * 17. + 22. + 32. + 0.5 * 37.
    assert math.isclose(total_flux, 81., rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(tarea, 3., rel_tol=1.e-8, abs_tol=1.e-8)

    image[:, 2] = np.nan

    (total_flux, f_var_poisson, f_var_rnoise, f_var_flat,
     bkg_flux, b_var_poisson, b_var_rnoise, b_var_flat,
     tarea, twht) = extract1d._extract_src_flux(
        image, var_poisson, var_rnoise, var_rflat,
        x, j, lam, srclim, weights, bkgmodels)

    assert np.isnan(total_flux)

    assert tarea == 0.


@pytest.mark.parametrize('test_type', ['all_empty', 'all_equal'])
def test_extract_src_flux_empty_interval(inputs_constant, test_type):
    (image, var_poisson, var_rnoise, var_rflat,
     x, j, lam, srclim, weights, bkgmodels) = inputs_constant

    if test_type == 'all_empty':
        # no limits provided
        srclim = []
    else:
        # empty extraction range: upper equals lower
        srclim[0][1] = srclim[0][0].copy()

    (total_flux, f_var_poisson, f_var_rnoise, f_var_flat,
     bkg_flux, b_var_poisson, b_var_rnoise, b_var_flat,
     tarea, twht) = extract1d._extract_src_flux(
        image, var_poisson, var_rnoise, var_rflat,
        x, j, lam, srclim, weights, bkgmodels)

    # empty interval, so no flux returned
    assert np.isnan(total_flux)
    assert bkg_flux == 0.
    assert tarea == 0.


@pytest.mark.parametrize('offset', [-100, 100])
def test_extract_src_flux_interval_out_of_range(inputs_constant, offset):
    (image, var_poisson, var_rnoise, var_rflat,
     x, j, lam, srclim, weights, bkgmodels) = inputs_constant

    # extraction limits out of range
    srclim[0][0] += offset
    srclim[0][1] += offset

    (total_flux, f_var_poisson, f_var_rnoise, f_var_flat,
     bkg_flux, b_var_poisson, b_var_rnoise, b_var_flat,
     tarea, twht) = extract1d._extract_src_flux(
        image, var_poisson, var_rnoise, var_rflat,
        x, j, lam, srclim, weights, bkgmodels)

    # empty interval, so no flux returned
    assert np.isnan(total_flux)
    assert bkg_flux == 0.
    assert tarea == 0.
