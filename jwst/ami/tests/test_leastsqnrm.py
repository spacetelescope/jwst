"""Unit tests for AMI leastsqnrm module."""

import numpy as np
import math
from numpy.testing import assert_allclose

from jwst.ami import leastsqnrm


def test_weighted_operations():
    img = np.arange(1, 5, dtype="f4").reshape((2, 2))
    model = np.arange(8, dtype="f4").reshape((2, 2, 2))
    dq = np.zeros((2, 2), dtype="int32")
    dq[1, 1] = 1
    x, res, cond, singvals = leastsqnrm.weighted_operations(img, model, dq)

    np.testing.assert_almost_equal(x, [-0.5, 1])
    assert cond is None
    np.testing.assert_almost_equal(singvals, [4.4870429, 0.1819676])


def test_leastsqnrm_rotatevectors():
    """Test of rotatevectors() in leastsqnrm module.
    Positive x decreases under slight rotation, and positive y
    increases under slight rotation.
    """
    vec = np.arange(8).reshape((4, 2)) + 1.0
    rot_vec = leastsqnrm.rotatevectors(vec, thetarad=0.001)

    true_rot_vec = np.array(
        [[0.9979995, 2.000999], [2.9959985, 4.002998], [4.9939975, 6.004997], [6.9919965, 8.006996]]
    )
    assert_allclose(rot_vec, true_rot_vec)

    rot_vec = leastsqnrm.rotatevectors(vec, thetarad=math.pi / 2.0)
    true_rot_vec = np.array([[-2.0, 1.0], [-4.0, 3.0], [-6.0, 5.0], [-8.0, 7.0]])
    assert_allclose(rot_vec, true_rot_vec)


def test_leastsqnrm_flip():
    """Test of flip() in leastsqnrm module.
    Change sign of 2nd coordinate of holes.
    """
    vec = np.arange(8).reshape((4, 2)) + 1.0
    flip_vec = leastsqnrm.flip(vec)

    true_flip_vec = np.array([[1.0, -2.0], [3.0, -4.0], [5.0, -6.0], [7.0, -8.0]])
    assert_allclose(flip_vec, true_flip_vec)


def test_leastsqnrm_mas2rad():
    """Test of mas2rad() in leastsqnrm module.
    Convert angle in milli arc-sec to radians.
    """
    mas = 1.0e8
    theta_rad = leastsqnrm.mas2rad(mas)

    true_theta_rad = mas * (10 ** (-3)) / (3600 * 180 / np.pi)
    assert_allclose(theta_rad, true_theta_rad)


def test_leastsqnrm_rad2mas():
    """Test of rad2mas() in leastsqnrm module.
    Convert input angle in radians to milli arc sec.
    """
    theta_rad = 1.0e-6
    mas = leastsqnrm.rad2mas(theta_rad)

    true_mas = theta_rad * (3600.0 * 180 / np.pi) * 10.0**3
    assert_allclose(mas, true_mas)


def test_leastsqnrm_sin2deltapistons():
    """Test of sin2deltapistons() in leastsqnrm module.
    Each baseline has one sine and one cosine fringe with a coefficient
    that depends on the piston difference between the two holes that make
    the baseline.  For a 7-hole mask there are 21 baselines and therefore
    there are 42 sine and cosine terms that contribute to the fringe model.
    This function calculates the sine of this piston difference.
    """
    # 4 holes (6 baselines) plus average flux per hole and DC offset
    coeffs = np.array([1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

    delta = leastsqnrm.sin2deltapistons(coeffs)

    true_delta = np.array([0.04849334, 0.08333333, 0.12340834])
    assert_allclose(delta, true_delta)


def test_leastsqnrm_cos2deltapistons():
    """Test of cos2deltapistons() in leastsqnrm module.
    Each baseline has one sine and one cosine fringe with a coefficient
    that depends on the piston difference between the two holes that make
    the baseline.  For a 7-hole mask there are 21 baselines and therefore
    there are 42 sine and cosine terms that contribute to the fringe model.
    This function calculate the cosine of this piston difference.
    """
    # 4 holes (6 baselines) plus average flux per hole and DC offset
    coeffs = np.array([1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

    delta = leastsqnrm.cos2deltapistons(coeffs)

    true_delta = np.array([0.21795289, 0.18450506, 0.14758362])
    assert_allclose(delta, true_delta)


def test_leastsqnrm_replacenan():
    """Test of replacenan() in leastsqnrm module.
    Replace singularities encountered in the analytical hexagon Fourier
    transform with the analytically derived limits. (pi/4)
    """
    arr = np.array([1.0, 5.6, np.nan, 5.3])

    rep_arr = leastsqnrm.replacenan(arr)

    true_rep_arr = np.array([1.0, 5.6, math.pi / 4.0, 5.3])
    assert_allclose(rep_arr, true_rep_arr)


def test_leastsqnrm_ffc():
    """Test of ffc in leastsqnrm module.
    Calculate cosine terms of analytic model.
    """
    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    ky = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    vv = np.arange(ASIZE)

    for ii in np.arange(ASIZE):
        kx[:, ii] = vv
        ky[ii, :] = vv

    # Clean up any attributes that may have been added earlier
    ffc = leastsqnrm.ffc
    for kk in list((ffc.__dict__).keys()):
        delattr(ffc, kk)

    ffc.N = 7
    ffc.lam = 2.3965000082171173e-06
    ffc.offx = 28.0
    ffc.offy = 28.0
    ffc.over = 3
    ffc.pitch = 9.099800633275124e-08
    ffc.ri = np.array([-0.01540951, -2.63995503])
    ffc.rj = np.array([-2.28627105, 0.01334504])
    ffc.size = (57, 57)

    ffc_arr = ffc(kx, ky)

    true_ffc_arr = np.array(
        [
            [-1.66542264, -0.68760215, 0.55667544, 1.58523217],
            [-1.99797288, -1.55759213, -0.51361892, 0.72939004],
            [-1.75826723, -1.98145931, -1.43680344, -0.3353627],
            [-1.01496177, -1.83780038, -1.94846124, -1.30406145],
        ]
    )

    assert_allclose(ffc_arr, true_ffc_arr)

    for kk in list((ffc.__dict__).keys()):
        delattr(ffc, kk)


def test_leastsqnrm_ffs():
    """Test of ffs in leastsqnrm module.
    Calculate sine terms of analytic model.
    """
    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    ky = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    vv = np.arange(ASIZE)

    for ii in np.arange(ASIZE):
        kx[:, ii] = vv
        ky[ii, :] = vv

    ffs = leastsqnrm.ffs
    for kk in list((ffs.__dict__).keys()):
        delattr(ffs, kk)

    ffs.N = 7
    ffs.lam = 2.3965000082171173e-06
    ffs.offx = 28.0
    ffs.offy = 28.0
    ffs.over = 3
    ffs.pitch = 9.099800633275124e-08
    ffs.ri = np.array([-0.01540951, -2.63995503])
    ffs.rj = np.array([-2.28627105, 0.01334504])
    ffs.size = (57, 57)

    ffs_arr = ffs(kx, ky)

    true_ffs_arr = np.array(
        [
            [-1.10741476, -1.878085, -1.92096654, -1.21944207],
            [-0.09002431, -1.2545544, -1.93292411, -1.86225406],
            [0.95315074, -0.27169652, -1.39125694, -1.97168249],
            [1.72332603, 0.7889802, -0.45110839, -1.51638509],
        ]
    )

    assert_allclose(ffs_arr, true_ffs_arr)

    for kk in list((ffs.__dict__).keys()):
        delattr(ffs, kk)


def test_leastsqnrm_closure_amplitudes():
    """Test of closure_amplitudes() in leastsqnrm module.
    Calculate the closure amplitudes.
    """
    amps = np.array([0.1, 0.2, 0.3, 1.0, 0.9, 0.5, 1.1, 0.7, 0.1, 1.0])

    n = 5  # number of holes

    cas = leastsqnrm.closure_amplitudes(amps, n=n)

    true_cas = np.array([0.7, 0.04545455, 0.3030303, 6.66666667, 18.0])
    assert_allclose(cas, true_cas, atol=1e-7)


def test_leastsqnrm_closurephase():
    """Test of closurephase in leastsqnrm module.
    Calculate closure phases between each pair of holes.
    """
    n = 7  # number of holes
    deltap = np.array(
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

    cps = leastsqnrm.closurephase(deltap, n=n)

    true_cps = np.array(
        [0.25, 0.5, 0.0, 0.07, 0.3, -0.8, 0.55, 0.03, 0.75, -0.35, -0.35, -0.65, 0.8, 1.2, 0.0]
    )
    assert_allclose(cps, true_cps, atol=1e-8)


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
