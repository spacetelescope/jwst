"""Unit tests for AMI leastsqnrm module."""

import pytest
import numpy as np
import math
from numpy.testing import assert_allclose

from jwst.ami import leastsqnrm


@pytest.mark.parametrize("pass_dq", [False, True])
def test_weighted_operations(pass_dq):
    img = np.arange(1, 5, dtype="f4").reshape((2, 2))
    model = np.arange(8, dtype="f4").reshape((2, 2, 2))
    if pass_dq:
        dq = np.zeros((2, 2), dtype="int32")
        dq[1, 1] = 1
    else:
        dq = None
    x, res, cond, singvals = leastsqnrm.weighted_operations(img, model, dq)

    np.testing.assert_almost_equal(x, [-0.5, 1])
    assert cond is None


@pytest.mark.parametrize("pass_dq", [False, True])
def test_matrix_operations(pass_dq):
    img = np.arange(1, 5, dtype="f4").reshape((2, 2))
    model = np.arange(8, dtype="f4").reshape((2, 2, 2))
    if pass_dq:
        dq = np.zeros((2, 2), dtype="int32")
        dq[1, 1] = 1
    else:
        dq = None
    x, res, cond, linfit_result = leastsqnrm.matrix_operations(img, model, dqm=dq, linfit=True)
    np.testing.assert_almost_equal(x, [-0.5, 1])
    assert isinstance(cond, float)


def test_leastsqnrm_replacenan():
    """Test of replacenan() in leastsqnrm module.
    Replace singularities encountered in the analytical hexagon Fourier
    transform with the analytically derived limits. (pi/4)
    """
    arr = np.array([1.0, 5.6, np.nan, 5.3])

    rep_arr = leastsqnrm.replacenan(arr)

    true_rep_arr = np.array([1.0, 5.6, math.pi / 4.0, 5.3])
    assert_allclose(rep_arr, true_rep_arr)


def test_leastsqnrm_closure_amplitudes():
    """Test of closure_amplitudes() in leastsqnrm module.
    Calculate the closure amplitudes.
    """
    amps = np.array([0.1, 0.2, 0.3, 1.0, 0.9, 0.5, 1.1, 0.7, 0.1, 1.0])

    n = 5  # number of holes

    cas = leastsqnrm.closure_amplitudes(amps, n=n)

    true_cas = np.array([0.7, 0.04545455, 0.3030303, 6.66666667, 18.0])
    assert_allclose(cas, true_cas, atol=1e-7)


def test_leastsqnrm_redundant_cps():
    """Test of redundant_cps in leastsqnrm module.
    Calculate closure phases for each set of 3 holes.
    """
    n = 7  # number of holes
    deltaps = np.array(
        [
            0.1,
            -0.2,
            0.3,
            0.2,
            0.05,
            -0.7,
            -0.05,
            0.7,
            0.1,
            0.02,
            -0.5,
            -0.05,
            0.7,
            0.1,
            0.3,
            0.4,
            -0.2,
            -0.3,
            0.2,
            0.5,
            0.3,
        ]
    )

    cps = leastsqnrm.redundant_cps(deltaps, n=n)

    true_cps = np.array(
        [
            0.25,
            0.5,
            0.0,
            0.07,
            0.3,
            -0.55,
            0.3,
            -0.15,
            0.8,
            0.5,
            0.05,
            0.7,
            0.35,
            1.4,
            1.05,
            -0.8,
            0.55,
            0.03,
            0.75,
            1.0,
            0.48,
            0.9,
            0.28,
            1.1,
            0.82,
            -0.35,
            -0.35,
            -0.65,
            0.8,
            0.9,
            0.1,
            0.8,
            1.2,
            0.4,
            0.0,
        ]
    )
    assert_allclose(cps, true_cps, atol=1e-8)


def test_leastsqnrm_populate_symmamparray():
    """Test of populate_symmamparray in leastsqnrm module.
    Populate the symmetric fringe amplitude array.
    """
    amps = np.array([0.1, 0.2, 0.3, 0.2, 0.05, 0.7, 0.3, 0.1, 0.2, 0.8])
    n = 5

    arr = leastsqnrm.populate_symmamparray(amps, n=n)

    true_arr = np.array(
        [
            [0.0, 0.1, 0.2, 0.3, 0.2],
            [0.1, 0.0, 0.05, 0.7, 0.3],
            [0.2, 0.05, 0.0, 0.1, 0.2],
            [0.3, 0.7, 0.1, 0.0, 0.8],
            [0.2, 0.3, 0.2, 0.8, 0.0],
        ]
    )
    assert_allclose(arr, true_arr, atol=1e-8)


def test_leastsqnrm_populate_antisymmphasearray():
    """Test of populate_antisymmphasearray in leastsqnrm module.
    Populate the antisymmetric fringe phase array.
    """
    deltaps = np.array([0.1, 0.2, 0.3, 0.2, 0.05, 0.7, 0.3, 0.1, 0.2, 0.8])
    n = 5

    arr = leastsqnrm.populate_antisymmphasearray(deltaps, n=n)

    true_arr = np.array(
        [
            [0.0, 0.1, 0.2, 0.3, 0.2],
            [-0.1, 0.0, 0.05, 0.7, 0.3],
            [-0.2, -0.05, 0.0, 0.1, 0.2],
            [-0.3, -0.7, -0.1, 0.0, 0.8],
            [-0.2, -0.3, -0.2, -0.8, 0.0],
        ]
    )
    assert_allclose(arr, true_arr, atol=1e-8)


def test_leastsqnrm_tan2visibilities():
    """Test of tan2visibilities in leastsqnrm module.
    From the solution to the fit, calculate the fringe amplitude and phase.
    """
    test_res = []  # to accumulate subtest comparisons

    coeffs = np.array([1.0, 0.2, -0.3, -0.1, 0.4, 0.2, -0.5, -0.1, 0.2, 0.4])

    amp, delta = leastsqnrm.tan2visibilities(coeffs)

    true_amp = np.array([0.36055513, 0.41231056, 0.53851648, 0.2236068])
    test_res.append(np.allclose(amp, true_amp, rtol=1.0e-7))

    true_delta = np.array([-0.98279372, 1.81577499, -1.19028995, 2.03444394])
    test_res.append(np.allclose(delta, true_delta, rtol=1.0e-7))

    assert np.all(test_res)


def test_leastsqnrm_multiplyenv():
    """Test of multiplyenv in leastsqnrm module.
    Multiply the envelope by each fringe 'image'.
    """
    env = np.array(
        [
            [1.00e-06, 1.00e-06, 1.01e-04],
            [1.00e-06, 1.00e-06, 1.01e-04],
            [1.10e-05, 1.10e-05, 1.00e-06],
        ]
    )

    ft = np.array(
        [
            [[4.0, 4.0, 4.0], [4.0, 4.0, 4.0], [4.0, 4.0, 4.0]],
            [[0.3, 0.3, 0.3], [0.3, 0.3, 0.4], [0.3, 0.3, 0.3]],
            [[-0.4, -0.4, -0.4], [-0.4, -0.4, -0.6], [-0.4, -0.4, -0.9]],
            [[0.8, 0.8, 0.8], [0.8, 0.8, 0.8], [1.0, 0.8, 1.3]],
        ]
    )

    # function is expecting a list, so make it one
    fringeterms = list(ft)

    full = leastsqnrm.multiplyenv(env, fringeterms)

    true_full = np.array(
        [
            [
                [4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e00],
                [4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e00],
                [4.04e-04, 3.03e-05, -4.04e-05, 8.08e-05, 1.00e00],
            ],
            [
                [4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e00],
                [4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e00],
                [4.04e-04, 4.04e-05, -6.06e-05, 8.08e-05, 1.00e00],
            ],
            [
                [4.40e-05, 3.30e-06, -4.40e-06, 1.10e-05, 1.00e00],
                [4.40e-05, 3.30e-06, -4.40e-06, 8.80e-06, 1.00e00],
                [4.00e-06, 3.00e-07, -9.00e-07, 1.30e-06, 1.00e00],
            ],
        ]
    )

    assert_allclose(full, true_full, atol=1e-8)
