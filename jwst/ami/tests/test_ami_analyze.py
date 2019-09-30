"""

Unit tests for ami_analyze

"""
import numpy as np

from jwst.ami import utils
from numpy.testing import assert_allclose

def test_utils_rebin():
    ''' test of rebin() and krebin() in utils module '''

    arr = np.arange(24).reshape((3,8))/10.
    rc = tuple((2,2))

    binned_arr = utils.rebin(arr, rc)

    true_arr = np.array([[5.1, 6.3, 7.5, 8.7]])
    assert_allclose(binned_arr, true_arr)


def test_utils_quadratic():
    ''' test of quadratic in utils module '''

    x = np.array( [0.5, 0.55, 0.55, 0.65, 0.70, 0.8, 0.85, 1., 1.01, 1.02,
                   1.03, 1.04, 1.05] )
    p = np.array( [-2., 3., 7.] )

    maxx, maxy, fit_vals = utils.quadratic( p, x )

    true_maxx = 0.75
    true_maxy = 8.125
    assert_allclose( [maxx, maxy], [true_maxx, true_maxy])

    true_fit_vals = np.array( [8., 8.045, 8.045, 8.105, 8.12, 8.12, 8.105, 8.,
                               7.9898, 7.9792, 7.9682, 7.9568, 7.945 ])
    assert_allclose( fit_vals, true_fit_vals )


def test_utils_findmax():
    ''' test of findmax in utils module '''

    mag = np.arange(9) +1.
    delt = 1.0E-7
    mag[2] += delt; mag[5] += delt  # Add a bit of deterministic noise
    mag[1] -= delt; mag[7] -= delt

    vals =( mag-3.)**2 +5   # Is quadratic ...
    vals[1] += delt; vals[6] += delt    # ... with a bit more noise
    vals[4] -= delt; vals[3] -= delt

    maxx, maxy = utils.findmax( mag, vals)

    true_maxx = 3.0
    true_maxy = 5.0
    assert_allclose( [maxx, maxy], [true_maxx, true_maxy])


def test_utils_makeA():
    ''' test of makeA in utils module '''

    nh = 4 # number of holes

    arr = utils.makeA( nh )

    true_arr = np.array([[-1.,  1.,  0.,  0.],
                         [-1.,  0.,  1.,  0.],
                         [-1.,  0.,  0.,  1.],
                         [ 0., -1.,  1.,  0.],
                         [ 0., -1.,  0.,  1.],
                         [ 0.,  0., -1.,  1.]])
    assert_allclose(arr, true_arr)


def test_utils_fringes2pistons():
    ''' test of fringes2pistons in utils module  '''

    fringephases = np.array([ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                              0.08, 0.09, 0.1 ])
    nholes = 5

    result = utils.fringes2pistons( fringephases, nholes )

    true_result = np.array([ 0.02, 0.034, 0.02, -0.014, -0.06 ])
    assert_allclose( result, true_result )


def test_utils_rcrosscorrelate():
    ''' test of rcrosscorrelate() in utils module  '''

    a = np.array([[ 2.,  3.,  4., 5.],
                  [ 6.,  7.,  8., 9.],
                  [10., 11., 12., 13.],
                  [14., 15., 16., 17.]])

    b = np.array([[-5., -4., -3., -2.],
                  [-1.,  0.,  1.,  2.],
                  [ 3.,  4.,  5.,  6.],
                  [ 7.,  8.,  9., 10.]])

    result = utils.rcrosscorrelate( a, b )

    true_result = np.array([[0.19865015, 0.20767971, 0.23476836, 0.20767971],
                            [0.34312299, 0.35215254, 0.3792412 , 0.35215254],
                            [0.77654151, 0.78557106, 0.81265972, 0.78557106],
                            [0.34312299, 0.35215254, 0.3792412 , 0.35215254]])
    assert_allclose( result, true_result )


def test_utils_crosscorrelate():
    ''' test of crosscorrelate() in utils module '''

    a = np.array([[ 2.,  3.,  4., 5.],
                  [ 6.,  7.,  8., 9.],
                  [10., 11., 12., 13.],
                  [14., 15., 16., 17.]])

    b = np.array([[-5., -4., -3., -2.],
                  [-1.,  0.,  1.,  2.],
                  [ 3.,  4.,  5.,  6.],
                  [ 7.,  8.,  9., 10.]])

    result = utils.crosscorrelate( a, b )

    true_result = np.array([[176.+0.j, 184.+0.j, 208.+0.j, 184.+0.j],
                            [304.+0.j, 312.+0.j, 336.+0.j, 312.+0.j],
                            [688.+0.j, 696.+0.j, 720.+0.j, 696.+0.j],
                            [304.+0.j, 312.+0.j, 336.+0.j, 312.+0.j]])
    assert_allclose( result, true_result )
