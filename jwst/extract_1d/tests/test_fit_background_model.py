"""
Test for extract_1d._fit_background_model
"""
import math
from copy import deepcopy

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

    b_lower = np.zeros(shape[1], dtype=np.float64) + 3.5    # 4, inclusive
    b_upper = np.zeros(shape[1], dtype=np.float64) + 4.5    # 4, inclusive
    bkglim = [[b_lower, b_upper]]
    bkg_order = 0
    bkg_fit = "poly"

    return image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order


def test_fit_background_model(inputs_constant):

    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(*inputs_constant)

    assert math.isclose(bkg_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 2


def test_fit_background_mean(inputs_constant):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant
    bkg_fit = "mean"

    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 2


def test_fit_background_median(inputs_constant):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant
    bkg_fit = "median"

    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 2


def test_handles_nan(inputs_constant):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant
    image[:, 2] = np.nan

    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 0


def test_handles_one_value(inputs_constant):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant
    image[np.where(image == 22)] = np.nan   # During extraction, only two pixels are returned "normally"; set one to Nan

    # If only one data point is available, the polynomial fit is forced to 0
    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 0


@pytest.mark.parametrize('test_type', ['all_empty', 'all_equal'])
def test_handles_empty_interval(inputs_constant, test_type):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant

    if test_type == 'all_empty':
        # no limits provided
        bkglim = []
    else:
        # empty extraction range: upper equals lower
        bkglim[0][1] = bkglim[0][0].copy()

    # No data available: background model is 0
    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 0


@pytest.mark.parametrize('offset', [-100, 100])
def test_handles_interval_out_of_range(inputs_constant, offset):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant
    bkglim[0][0] += offset
    bkglim[0][1] += offset

    # No data available: background model is 0
    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 0


def test_handles_one_empty_interval(inputs_constant):
    image, var_poisson, var_rnoise, var_rflat, x, j, bkglim, bkg_fit, bkg_order = inputs_constant

    # add an extra interval that is empty
    bkglim.append(deepcopy(bkglim[0]))
    bkglim[1][0] += 2.0
    bkglim[1][1] = bkglim[1][0].copy()
    print(bkglim)

    # should ignore the second interval and return a valid answer for the first
    (bkg_model, b_var_poisson_model, b_var_rnoise_model, b_var_flat_model, npts) = \
        extract1d._fit_background_model(image, var_poisson, var_rnoise, var_rflat,
                                        x, j, bkglim, bkg_fit, bkg_order)

    assert math.isclose(bkg_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_poisson_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_poisson_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_rnoise_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_rnoise_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert math.isclose(b_var_flat_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(b_var_flat_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 2
