"""Unit tests for AMI utils module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
import stdatamodels.jwst.datamodels as dm

from jwst.ami import utils


@pytest.mark.parametrize(
    "shape, center",
    [
        ((10, 10), (4.5, 4.5)),
        ((11, 11), (5, 5)),
    ],
)
def test_centerpoint(shape, center):
    assert utils.centerpoint(shape) == center


def test_find_centroid():
    arr = np.zeros((30, 30), dtype="f4")
    arr[15, 15] = 1
    assert_allclose(utils.find_centroid(arr), (0.5, 0.5))


@pytest.mark.parametrize(
    "mas, rad",
    [
        (206264.8062471, 0.001),
        (103132403.12355, 0.5),
    ],
)
def test_mas2rad(mas, rad):
    assert np.isclose(utils.mas2rad(mas), rad)


def test_utils_rebin():
    """Test of rebin() and krebin() in utils module"""
    arr = np.arange(24).reshape((3, 8)) / 10.0
    rc = tuple((2, 2))

    binned_arr = utils.rebin(arr, rc)

    true_arr = np.array([[5.1, 6.3, 7.5, 8.7]])
    assert_allclose(binned_arr, true_arr)


def test_utils_quadratic():
    """Test of quadratic in utils module"""
    x = np.array([0.5, 0.55, 0.55, 0.65, 0.70, 0.8, 0.85, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05])
    p = np.array([-2.0, 3.0, 7.0])

    maxx, maxy, fit_vals = utils.quadratic(p, x)

    true_maxx = 0.75
    true_maxy = 8.125
    assert_allclose([maxx, maxy], [true_maxx, true_maxy])

    true_fit_vals = np.array(
        [8.0, 8.045, 8.045, 8.105, 8.12, 8.12, 8.105, 8.0, 7.9898, 7.9792, 7.9682, 7.9568, 7.945]
    )
    assert_allclose(fit_vals, true_fit_vals)


def test_utils_findmax():
    """Test of findmax in utils module"""
    mag = np.arange(9) + 1.0
    delt = 1.0e-7
    mag[2] += delt
    mag[5] += delt  # Add a bit of deterministic noise
    mag[1] -= delt
    mag[7] -= delt

    vals = (mag - 3.0) ** 2 + 5  # Is quadratic ...
    vals[1] += delt
    vals[6] += delt  # ... with a bit more noise
    vals[4] -= delt
    vals[3] -= delt

    maxx, maxy = utils.findmax(mag, vals)

    true_maxx = 3.0
    true_maxy = 5.0
    assert_allclose([maxx, maxy], [true_maxx, true_maxy])


def test_utils_makeA():
    """Test of makeA in utils module"""
    nh = 4  # number of holes
    arr = utils.make_a(nh)

    true_arr = np.array(
        [
            [-1.0, 1.0, 0.0, 0.0],
            [-1.0, 0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0, 1.0],
            [0.0, -1.0, 1.0, 0.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, -1.0, 1.0],
        ]
    )
    assert_allclose(arr, true_arr)


def test_utils_fringes2pistons():
    """Test of fringes2pistons in utils module"""
    fringephases = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    nholes = 5

    result = utils.fringes2pistons(fringephases, nholes)

    true_result = np.array([-0.02, -0.034, -0.02, 0.014, 0.06])
    assert_allclose(result, true_result)


def test_utils_rcrosscorrelate():
    """Test of rcrosscorrelate() in utils module"""
    a = np.array(
        [
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0, 13.0],
            [14.0, 15.0, 16.0, 17.0],
        ]
    )

    b = np.array(
        [
            [-5.0, -4.0, -3.0, -2.0],
            [-1.0, 0.0, 1.0, 2.0],
            [3.0, 4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0, 10.0],
        ]
    )

    result = utils.rcrosscorrelate(a, b)

    true_result = np.array(
        [
            [0.19865015, 0.20767971, 0.23476836, 0.20767971],
            [0.34312299, 0.35215254, 0.3792412, 0.35215254],
            [0.77654151, 0.78557106, 0.81265972, 0.78557106],
            [0.34312299, 0.35215254, 0.3792412, 0.35215254],
        ]
    )
    assert_allclose(result, true_result)


def test_utils_crosscorrelate():
    """Test of crosscorrelate() in utils module"""
    a = np.array(
        [
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0, 13.0],
            [14.0, 15.0, 16.0, 17.0],
        ]
    )

    b = np.array(
        [
            [-5.0, -4.0, -3.0, -2.0],
            [-1.0, 0.0, 1.0, 2.0],
            [3.0, 4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0, 10.0],
        ]
    )

    result = utils.crosscorrelate(a, b)

    true_result = np.array(
        [
            [176.0 + 0.0j, 184.0 + 0.0j, 208.0 + 0.0j, 184.0 + 0.0j],
            [304.0 + 0.0j, 312.0 + 0.0j, 336.0 + 0.0j, 312.0 + 0.0j],
            [688.0 + 0.0j, 696.0 + 0.0j, 720.0 + 0.0j, 696.0 + 0.0j],
            [304.0 + 0.0j, 312.0 + 0.0j, 336.0 + 0.0j, 312.0 + 0.0j],
        ]
    )
    assert_allclose(result, true_result)


def test_utils_imgmedian():
    """Test of img_median_replace() in utils module"""
    # create input image model containing NaN's and DO_NOT_USE flags
    data = np.array(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [6.0, 0.0, 8.0, 9.0, 10.0],
            [11.0, 12.0, 13.0, 14.0, 15.0],
            [16.0, 17.0, np.nan, 0.0, 20.0],
            [21.0, 22.0, 23.0, 24.0, 25.0],
        ],
        dtype=np.float32,
    )

    dq = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 1, 4, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0],
        ],
        dtype=np.uint32,
    )

    input_model = dm.ImageModel(data=data, dq=dq)

    # send to img_median_replace to replace bad pixels
    input_model = utils.img_median_replace(input_model, box_size=3)

    expected_result = np.array(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0, 10.0],
            [11.0, 12.0, 13.0, 14.0, 15.0],
            [16.0, 17.0, 17.0, 18.5, 20.0],
            [21.0, 22.0, 23.0, 24.0, 25.0],
        ],
        dtype=np.float32,
    )

    assert_allclose(input_model.data, expected_result)
