"""

Unit tests for ami_analyze

"""
import numpy as np
import math

from jwst import datamodels
from jwst.ami import utils, leastsqnrm, hexee, webb_psf
from jwst.ami.leastsqnrm import hexpb, ffc, ffs, return_CAs
from jwst.ami.leastsqnrm import closurephase, redundant_cps
from jwst.ami.leastsqnrm import populate_symmamparray
from jwst.ami.leastsqnrm import populate_antisymmphasearray
from jwst.ami.leastsqnrm import tan2visibilities, model_array
from jwst.ami.analyticnrm2 import interf, PSF, phasor, ASFhex

from numpy.testing import assert_allclose

#---------------------------------------------------------------
# utils module tests:
#
def test_utils_rebin():
    ''' Test of rebin() and krebin() in utils module '''

    arr = np.arange(24).reshape((3,8))/10.
    rc = tuple((2,2))

    binned_arr = utils.rebin(arr, rc)

    true_arr = np.array([[5.1, 6.3, 7.5, 8.7]])
    assert_allclose(binned_arr, true_arr)


def test_utils_quadratic():
    ''' Test of quadratic in utils module '''

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
    ''' Test of findmax in utils module '''

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
    ''' Test of makeA in utils module '''

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
    ''' Test of fringes2pistons in utils module  '''

    fringephases = np.array([ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                              0.08, 0.09, 0.1 ])
    nholes = 5

    result = utils.fringes2pistons( fringephases, nholes )

    true_result = np.array([ 0.02, 0.034, 0.02, -0.014, -0.06 ])
    assert_allclose( result, true_result )


def test_utils_rcrosscorrelate():
    ''' Test of rcrosscorrelate() in utils module '''

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
    ''' Test of crosscorrelate() in utils module '''

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


#---------------------------------------------------------------
# leastsqnrm module tests:
#
def test_leastsqnrm_rotatevectors():
    ''' Test of rotatevectors() in leastsqnrm module.
        Positive x decreases under slight rotation, and positive y
        increases under slight rotation.
    '''
    vec = np.arange(8).reshape((4,2)) + 1.
    rot_vec = leastsqnrm.rotatevectors(vec, thetarad=0.001)

    true_rot_vec = np.array([[0.9979995, 2.000999],
                             [2.9959985, 4.002998],
                             [4.9939975, 6.004997],
                             [6.9919965, 8.006996]])
    assert_allclose( rot_vec, true_rot_vec )

    rot_vec = leastsqnrm.rotatevectors(vec, thetarad=math.pi/2.)
    true_rot_vec = np.array([[-2., 1.],
                             [-4., 3.],
                             [-6., 5.],
                             [-8., 7.]])
    assert_allclose( rot_vec, true_rot_vec )


def test_leastsqnrm_flip():
    ''' Test of flip() in leastsqnrm module.
        Change sign of 2nd coordinate of holes.
    '''
    vec = np.arange(8).reshape((4,2)) + 1.
    flip_vec = leastsqnrm.flip(vec)

    true_flip_vec = np.array([[ 1., -2.],
                              [ 3., -4.],
                              [ 5., -6.],
                              [ 7., -8.]])
    assert_allclose( flip_vec, true_flip_vec )


def test_leastsqnrm_mas2rad():
    ''' Test of mas2rad() in leastsqnrm module.
        Convert angle in milli arc-sec to radians.
    '''
    mas = 1.E8
    theta_rad = leastsqnrm.mas2rad( mas )

    true_theta_rad = mas * (10**(-3)) / (3600 * 180 / np.pi)
    assert_allclose( theta_rad, true_theta_rad )


def test_leastsqnrm_rad2mas():
    ''' Test of rad2mas() in leastsqnrm module.
        Convert input angle in radians to milli arc sec.
    '''
    theta_rad = 1.E-6
    mas = leastsqnrm.rad2mas( theta_rad )

    true_mas = theta_rad * (3600. * 180 / np.pi) * 10.**3
    assert_allclose( mas, true_mas )


def test_leastsqnrm_sin2deltapistons():
    ''' Test of sin2deltapistons() in leastsqnrm module.
        Each baseline has one sine and one cosine fringe with a coefficient
        that depends on the piston difference between the two holes that make
        the baseline.  For a 7-hole mask there are 21 baselines and therefore
        there are 42 sine and cosine terms that contribute to the fringe model.
        This function calculates the sine of this piston difference.
    '''
    # 4 holes (6 baselines) plus average flux per hole and DC offset
    coeffs = np.array([1., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

    delta = leastsqnrm.sin2deltapistons( coeffs )

    true_delta = np.array([0.04849334, 0.08333333, 0.12340834])
    assert_allclose( delta, true_delta)


def test_leastsqnrm_cos2deltapistons():
    ''' Test of cos2deltapistons() in leastsqnrm module.
        Each baseline has one sine and one cosine fringe with a coefficient
        that depends on the piston difference between the two holes that make
        the baseline.  For a 7-hole mask there are 21 baselines and therefore
        there are 42 sine and cosine terms that contribute to the fringe model.
        This function calculate the cosine of this piston difference.
    '''
    # 4 holes (6 baselines) plus average flux per hole and DC offset
    coeffs = np.array([1., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

    delta = leastsqnrm.cos2deltapistons( coeffs )

    true_delta = np.array([0.21795289, 0.18450506, 0.14758362])
    assert_allclose( delta, true_delta)


def test_leastsqnrm_replacenan():
    ''' Test of replacenan() in leastsqnrm module.
        Replace singularities encountered in the analytical hexagon Fourier
        transform with the analytically derived limits. (pi/4)
    '''
    arr = np.array([ 1., 5.6, np.nan, 5.3 ])

    rep_arr = leastsqnrm.replacenan( arr )

    true_rep_arr = np.array([1., 5.6, math.pi/4., 5.3 ])
    assert_allclose( rep_arr, true_rep_arr)


def test_leastsqnrm_hexpb():
    ''' Test of hexpb() in leastsqnrm module.
        Calculate the primary beam for hexagonal holes.
    '''
    # Clean up any attributes that may have been added earlier
    for kk in list( (hexpb.__dict__).keys()):
        delattr( hexpb, kk )

    hexpb.d = 0.5
    hexpb.lam = 2.0e-06
    hexpb.offx = 28.0
    hexpb.offy = 28.0
    hexpb.pitch = 1.0e-07
    hexpb.shape = 'hex'
    hexpb.size =(3, 3)

    hexpb_arr = hexpb()

    true_hexpb_arr = np.array( [[0.01520087, 0.01901502, 0.02328432],
                                [0.01912038, 0.02356723, 0.02850747],
                                [0.02349951, 0.02861771, 0.03426836]] )
    assert_allclose( hexpb_arr, true_hexpb_arr, atol=1E-7)

    # Clean up attributes that have been added
    for kk in list( (hexpb.__dict__).keys()):
        delattr( hexpb, kk )


def test_leastsqnrm_model_array():
    ''' Test of model_array in leastsqnrm module.
        Create a model using the specified wavelength.
    '''
    import warnings

    test_res = [] # to accumulate subtest comparisons

    modelctrs = np.array([[-0.01540951, -2.63995503],
                          [-2.28627105,  0.01334504],
                          [ 2.2785663,  -1.33332266],
                          [-2.2785663,   1.33332266],
                          [-1.1315734,   1.98663876],
                          [ 2.29397581,  1.30663257],
                          [ 1.15468766,  1.97329378]])
    lam = 2.3965000082171173e-06
    oversample = 3  # oversample factor
    modelpix = 2.729940189982537e-07
    fov = 19
    hole_d = 0.8 # hole diameter

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "invalid value*", RuntimeWarning)
        warnings.filterwarnings("ignore", "divide by zero*", RuntimeWarning)
        pb, ff = model_array(modelctrs, lam, oversample, modelpix, fov,
                              hole_d, centering="PIXELCENTERED", shape="hex")

    pb_1 = pb[15:18,20:23]
    true_pb_1 = np.array([[0.40298655, 0.42031716, 0.43582878],
                          [0.43220726, 0.45055339, 0.46696826],
                          [0.46060077, 0.47992491, 0.49720953]])
    test_res.append( np.allclose( pb_1, true_pb_1, rtol=1.e-7 ))

    pb_2 = pb[25:28,10:13]
    true_pb_2 = np.array([[0.30173028, 0.33450152, 0.36798482],
                          [0.30622563, 0.33939216, 0.37327326],
                          [0.30894849, 0.34235409, 0.37647573]])
    test_res.append( np.allclose( pb_2, true_pb_2, rtol=1.e-7 ))

    ff_1 = ff[1][15:18,20:23] # slice 1
    true_ff_1 = np.array([[-1.45863979, -0.5441426,   0.52620696],
                          [-1.98551072, -1.57725007, -0.71723623],
                          [-1.74296371, -1.9991473,  -1.68273863]])
    test_res.append( np.allclose( ff_1, true_ff_1, rtol=1.e-7 ))

    ff_5 = ff[5][25:28,10:13] # slice 5
    true_ff_5 = np.array([[ 1.65967899,  1.99729328, 1.76662698],
                          [ 0.06173976,  1.0806429,  1.79207573],
                          [-1.58764681, -0.73650014, 0.32419949]])
    test_res.append( np.allclose( ff_5, true_ff_5, rtol=1.e-7 ))

    assert np.alltrue( test_res )


def test_leastsqnrm_ffc():
    ''' Test of ffc in leastsqnrm module.
        Calculate cosine terms of analytic model.
    '''
    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    ky = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    vv = np.arange(ASIZE)

    for ii in np.arange(ASIZE):
        kx[:,ii] = vv
        ky[ii,:] = vv

    # Clean up any attributes that may have been added earlier
    for kk in list( (ffc.__dict__).keys()):
        delattr( ffc, kk )

    ffc.N = 7
    ffc.lam = 2.3965000082171173e-06
    ffc.offx = 28.0
    ffc.offy = 28.0
    ffc.over = 3
    ffc.pitch = 9.099800633275124e-08
    ffc.ri = np.array([-0.01540951, -2.63995503])
    ffc.rj = np.array([-2.28627105,  0.01334504])
    ffc.size = (57, 57)

    ffc_arr = ffc( kx, ky)

    true_ffc_arr = np.array([[-1.66542264, -0.68760215,  0.55667544,  1.58523217],
                             [-1.99797288, -1.55759213, -0.51361892,  0.72939004],
                             [-1.75826723, -1.98145931, -1.43680344, -0.3353627 ],
                             [-1.01496177, -1.83780038, -1.94846124, -1.30406145]])

    assert_allclose( ffc_arr, true_ffc_arr )

    for kk in list( (ffc.__dict__).keys()):
        delattr( ffc, kk )


def test_leastsqnrm_ffs():
    ''' Test of ffs in leastsqnrm module.
        Calculate sine terms of analytic model.
    '''
    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE ).reshape((ASIZE, ASIZE))
    ky = np.arange(ASIZE * ASIZE ).reshape((ASIZE, ASIZE))
    vv = np.arange(ASIZE)

    for ii in np.arange(ASIZE):
        kx[:,ii] = vv
        ky[ii,:] = vv

    for kk in list( (ffs.__dict__).keys()):
        delattr( ffs, kk )

    ffs.N = 7
    ffs.lam = 2.3965000082171173e-06
    ffs.offx = 28.0
    ffs.offy = 28.0
    ffs.over = 3
    ffs.pitch = 9.099800633275124e-08
    ffs.ri = np.array([-0.01540951, -2.63995503])
    ffs.rj = np.array([-2.28627105,  0.01334504])
    ffs.size = (57, 57)

    ffs_arr = ffs( kx, ky)

    true_ffs_arr = np.array([[-1.10741476, -1.878085, -1.92096654, -1.21944207],
                             [-0.09002431, -1.2545544, -1.93292411, -1.86225406],
                             [ 0.95315074, -0.27169652, -1.39125694, -1.97168249],
                             [ 1.72332603, 0.7889802, -0.45110839, -1.51638509]])

    assert_allclose( ffs_arr, true_ffs_arr )

    for kk in list( (ffs.__dict__).keys()):
        delattr( ffs, kk )


def test_leastsqnrm_return_CAs():
    ''' Test of return_CAs in leastsqnrm module.
        Calculate the closure amplitudes.
    '''
    amps = np.array([ 0.1, 0.2, 0.3, 1.0, 0.9, 0.5 ,1.1, 0.7, 0.1, 1.0 ])

    N = 5 # number of holes

    CAs = return_CAs(amps, N=N)

    true_CAs = np.array([ 0.7, 0.04545455, 0.3030303, 6.66666667, 18.])
    assert_allclose( CAs, true_CAs, atol=1E-7 )


def test_leastsqnrm_closurephase():
    '''  Test of closurephase in leastsqnrm module.
         Calculate closure phases between each pair of holes.
    '''

    N = 7 # number of holes
    deltap = np.array([ 0.1, -0.2,  0.3, 0.2, 0.05,
                       -0.7, -0.05, 0.7, 0.1, 0.02,
                       -0.5, -0.05, 0.7, 0.1, 0.3,
                        0.4, -0.2, -0.3, 0.2, 0.5, 0.3 ])

    cps = closurephase( deltap, N=N )

    true_cps = np.array([ 0.25,  0.5,  0.0,  0.07 , 0.3,
                         -0.8,   0.55, 0.03, 0.75, -0.35,
                         -0.35, -0.65, 0.8,  1.2,   0.0 ])
    assert_allclose( cps, true_cps, atol=1E-8 )


def test_leastsqnrm_redundant_cps():
    ''' Test of redundant_cps in leastsqnrm module.
        Calculate closure phases for each set of 3 holes.
    '''
    N = 7 # number of holes
    deltaps = np.array([ 0.1, -0.2,  0.3, 0.2, 0.05,
                       -0.7, -0.05, 0.7, 0.1, 0.02,
                       -0.5, -0.05, 0.7, 0.1, 0.3,
                        0.4, -0.2, -0.3, 0.2, 0.5, 0.3 ])

    cps = redundant_cps( deltaps, N=N )

    true_cps = np.array([ 0.25,   0.5,   0.0,  0.07, 0.3,
                         -0.55,   0.3,  -0.15, 0.8,  0.5,
                          0.05,   0.7,   0.35, 1.4,  1.05,
                          -0.8,   0.55,  0.03, 0.75, 1.0,
                           0.48,  0.9,   0.28, 1.1,  0.82,
                          -0.35, -0.35, -0.65, 0.8,  0.9,
                           0.1,   0.8,   1.2,  0.4,  0.0 ])
    assert_allclose( cps, true_cps, atol=1E-8 )


def test_leastsqnrm_populate_symmamparray():
    ''' Test of populate_symmamparray in leastsqnrm module.
        Populate the symmetric fringe amplitude array.
    '''
    amps = np.array([ 0.1, 0.2,  0.3, 0.2, 0.05, 0.7, 0.3, 0.1, 0.2, 0.8 ])
    N = 5

    arr = populate_symmamparray( amps, N=N )

    true_arr = np.array([[0.0, 0.1,  0.2,  0.3,  0.2 ],
                         [0.1, 0.0,  0.05, 0.7,  0.3 ],
                         [0.2, 0.05, 0.0,  0.1,  0.2 ],
                         [0.3, 0.7,  0.1,  0.0,  0.8 ],
                         [0.2, 0.3,  0.2,  0.8,  0.0 ]])
    assert_allclose( arr, true_arr, atol=1E-8 )


def test_leastsqnrm_populate_antisymmphasearray():
    ''' Test of populate_antisymmphasearray in leastsqnrm module.
        Populate the antisymmetric fringe phase array.
    '''
    deltaps = np.array([ 0.1, 0.2,  0.3, 0.2, 0.05, 0.7, 0.3, 0.1, 0.2, 0.8 ])
    N = 5

    arr = populate_antisymmphasearray( deltaps, N=N )

    true_arr = np.array([[0.0, 0.1,  0.2,  0.3,  0.2 ],
                         [-0.1, 0.0,  0.05, 0.7,  0.3 ],
                         [-0.2, -0.05, 0.0,  0.1,  0.2 ],
                         [-0.3, -0.7,  -0.1,  0.0,  0.8 ],
                         [-0.2, -0.3,  -0.2,  -0.8,  0.0 ]])
    assert_allclose( arr, true_arr, atol=1E-8 )


def test_leastsqnrm_tan2visibilities():
    ''' Test of tan2visibilities in leastsqnrm module.
        From the solution to the fit, calculate the fringe amplitude and phase.
    '''
    test_res = [] # to accumulate subtest comparisons

    coeffs = np.array([ 1.0, 0.2, -0.3, -0.1, 0.4, 0.2, -0.5, -0.1, 0.2, 0.4])

    amp, delta = tan2visibilities(coeffs)

    true_amp = np.array([0.36055513, 0.41231056, 0.53851648, 0.2236068 ])
    test_res.append( np.allclose( amp, true_amp, rtol=1.e-7 ))

    true_delta = np.array([-0.98279372,  1.81577499, -1.19028995, 2.03444394])
    test_res.append( np.allclose( delta, true_delta, rtol=1.e-7 ))

    assert np.alltrue( test_res )


def test_leastsqnrm_multiplyenv():
    ''' Test of multiplyenv in leastsqnrm module.
        Multiply the envelope by each fringe 'image'.
    '''
    env = np.array([[1.00e-06, 1.00e-06, 1.01e-04],
                    [1.00e-06, 1.00e-06, 1.01e-04],
                    [1.10e-05, 1.10e-05, 1.00e-06]])

    ft = np.array([[[ 4., 4., 4. ],
                    [ 4., 4., 4. ],
                    [ 4., 4., 4. ]],
                   [[ 0.3, 0.3, 0.3],
                    [ 0.3, 0.3, 0.4],
                    [ 0.3, 0.3, 0.3]],
                   [[-0.4, -0.4, -0.4],
                    [-0.4, -0.4, -0.6],
                    [-0.4, -0.4, -0.9]],
                   [[ 0.8, 0.8, 0.8],
                    [ 0.8, 0.8, 0.8],
                    [ 1., 0.8, 1.3]]])

   # function is expecting a list, so make it one
    fringeterms = list( ft )

    full = leastsqnrm.multiplyenv( env, fringeterms )

    true_full = np.array([[[ 4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e+00],
                           [ 4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e+00],
                           [ 4.04e-04, 3.03e-05, -4.04e-05, 8.08e-05, 1.00e+00]],
                          [[ 4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e+00],
                           [ 4.00e-06, 3.00e-07, -4.00e-07, 8.00e-07, 1.00e+00],
                           [ 4.04e-04, 4.04e-05, -6.06e-05, 8.08e-05, 1.00e+00]],
                          [[ 4.40e-05, 3.30e-06, -4.40e-06, 1.10e-05, 1.00e+00],
                           [ 4.40e-05, 3.30e-06, -4.40e-06, 8.80e-06, 1.00e+00],
                           [ 4.00e-06, 3.00e-07, -9.00e-07, 1.30e-06, 1.00e+00]]])

    assert_allclose( full, true_full, atol=1E-8 )


#---------------------------------------------------------------
# hexee module tests:
#
def test_hexee_g_eeAG():
    ''' Test of g_eeAG() in the hexee module.
        Calculate the Fourier transform of one half of a hexagon that is
        bisected from one corner to its diametrically opposite corner.
    '''
    xi, eta, kwargs = setup_hexee()
    g = hexee.g_eeAG(xi, eta, **kwargs)

    true_g = np.array([[ -0.04454286+0.05015766j, -0.04164985+0.06041733j,
                         -0.03830953+0.07099764j ],
                       [ -0.04072437+0.05375103j, -0.03729262+0.06415232j,
                         -0.03340318+0.07486623j ],
                       [ -0.03657856+0.05703437j, -0.03258885+0.06754246j,
                         -0.02813134+0.07835476j ]])

    assert_allclose(g, true_g)


def test_hexee_glimit():
    ''' Test of glimit() in the hexee module.
        Calculate the analytic limit of the Fourier transform of one half of the
        hexagon along eta=0.
    '''
    xi, eta, kwargs = setup_hexee()
    g = hexee.glimit(xi, eta, **kwargs)

    true_g = np.array( [[0.07105571+0.28088478j, 0.07105571+0.28088478j,
                         0.07105571+0.28088478j ],
                        [0.08609692+0.28598645j, 0.08609692+0.28598645j,
                         0.08609692+0.28598645j ],
                        [0.10178022+0.29008864j, 0.10178022+0.29008864j,
                         0.10178022+0.29008864j]] )

    assert_allclose( g, true_g )


#---------------------------------------------------------------
# analyticnrm2 module tests:
#
def test_analyticnrm2_PSF():
    ''' Test of PSF() in the analyticnrm2 module '''

    pixel, fov, oversample, ctrs, d, lam, phi, centering = setup_SF()

    shape = "hex"

    psf =  PSF( pixel, fov, oversample, ctrs, d, lam, phi, centering = centering, shape = shape)

    true_psf = np.array(
        [[ 0.55676539, 5.13790691, 9.32707852, 7.95659204, 3.05720862, 0.87093812],
         [ 2.883065, 13.4868413, 22.44612867, 20.24062814, 9.2865436, 0.8260206 ],
         [ 4.66981271, 19.68183766, 33.25399012, 32.22248535, 17.64512595, 3.49243422],
         [ 3.49243422, 17.64512595, 32.22248535, 33.25399012, 19.68183766, 4.66981271],
         [ 0.8260206, 9.2865436, 20.24062814, 22.44612867, 13.4868413, 2.883065 ],
         [ 0.87093812, 3.05720862, 7.95659204, 9.32707852, 5.13790691, 0.55676539]]
         )

    assert_allclose(psf, true_psf, atol=1E-7 )


def test_analyticnrm2_ASFhex():
    ''' Test of ASFhex() in the analyticnrm2 module FOR HEX '''

    pixel, fov, oversample, ctrs, d, lam, phi, centering = setup_SF()

    asf = ASFhex(pixel, fov, oversample, ctrs, d, lam, phi, centering)

    true_asf = np.array(
        [[ 0.74611461-0.00885311j, 2.14239062+0.74031707j, 2.7532696 +1.32158428j,
           2.35075667+1.55901736j, 1.09507456+1.3630922j, -0.53224117+0.7665882j ],
         [ 1.69794671-0.00647795j, 3.6397516 +0.48892695j, 4.66263759+0.840202j,
           4.39936028+0.9414124j,  2.95391668+0.74894583j, 0.85814396+0.29934852j],
         [ 2.16095243+0.00986439j, 4.43309812+0.17169358j, 5.75956679+0.28527278j,
           5.66797969+0.31063092j, 4.19436059+0.22905275j, 1.86810611+0.05112528j],
         [ 1.86810611-0.05112528j, 4.19436059-0.22905275j, 5.66797969-0.31063092j,
           5.75956679-0.28527278j, 4.43309812-0.17169358j, 2.16095243-0.00986439j],
         [ 0.85814396-0.29934852j, 2.95391668-0.74894583j, 4.39936028-0.9414124j,
           4.66263759-0.840202j,   3.6397516 -0.48892695j, 1.69794671+0.00647795j],
         [-0.53224117-0.7665882j,  1.09507456-1.3630922j,  2.35075667-1.55901736j,
           2.7532696 -1.32158428j, 2.14239062-0.74031707j, 0.74611461+0.00885311j]]
    )

    assert_allclose(asf, true_asf, atol=1E-7 )


def test_analyticnrm2_interf():
    ''' Test of interf() in the analyticnrm2 module '''

    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    ky = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    vv = np.arange(ASIZE)

    for ii in np.arange(ASIZE):
        kx[:,ii] = vv
        ky[ii,:] = vv

    # Clean up any attributes that may have been added earlier
    for kk in list( (interf.__dict__).keys()):
        delattr( interf, kk )

    pixel, fov, oversample, ctrs, d, lam, phi, centering = setup_SF()

    interf.lam = lam
    interf.offx = 0.5
    interf.offy = 0.5
    interf.ctrs = ctrs
    interf.d = d
    interf.phi = phi
    interf.pitch = pixel / float(oversample)

    interference = interf(kx, ky)

    true_interference = np.array(
       [[6.65604548+0.32967559j, 6.55020282-0.35898074j,
         5.10088264-1.09153011j, 2.74363797-1.81957549j ],
        [6.55020282+0.35898074j, 6.65604548-0.32967559j,
         5.40614218-0.97418068j, 3.21342278-1.5424603j ],
        [4.86319369+0.26557752j, 5.14000033-0.19907186j,
         4.23407392-0.56876212j, 2.50871764-0.86690377j ],
        [2.18032216+0.05966983j, 2.52211181-0.01151302j,
         1.98827833+0.00758561j, 0.87949091+0.01043571j]]
        )

    assert_allclose(interference, true_interference, atol=1E-7 )


def test_analyticnrm2_phasor():
    ''' Test of phasor() in the analyticnrm2 module '''

    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))

    for ii in np.arange(ASIZE):
        kx[:,ii] = ii

    ky = kx.transpose()

    hx = 0.06864653345335156
    hy = -2.6391073592116028

    lam = 2.3965000082171173e-06
    phi = 0.0
    pitch = 1.0375012775744072e-07

    result = phasor( kx, ky, hx, hy, lam, phi, pitch )

    true_result = np.array(
      [[ 1. -0.j, 0.99982567-0.01867173j, 0.99930273-0.03733694j,
         0.99843138-0.05598914j ],
       [ 0.75320598+0.65778473j, 0.76535665+0.6436064j, 0.77724047+0.62920367j,
         0.78885329+0.61458156j ],
       [ 0.1346385+0.99089478j, 0.15311675+0.98820811j, 0.1715416+0.98517688j,
         0.18990665+0.98180215j ],
       [-0.55038493+0.83491103j, -0.53469975+0.84504211j, -0.51882814+0.85487856j,
        -0.50277564+0.86441695j ]]
      )

    assert_allclose( result, true_result, atol=1E-7 )


#---------------------------------------------------------------
# webb_psf module tests:
#
def test_webb_psf():
    ''' Test of PSF() in the webb_psf module:
        Load the filter bandpass data using files from WebbPSF
    '''
    throughput_reffile = "jwst_niriss_throughput_0005.fits"
    filter_model = datamodels.ThroughputModel( throughput_reffile )
    bindown = 80

    band = webb_psf.get_webbpsf_filter(filter_model, specbin=bindown)

    true_band = np.array( [[5.67496444e-01, 2.43628541e-06],
                           [9.63244237e-01, 2.52823216e-06],
                           [9.81308160e-01, 2.62091254e-06],
                           [9.91945652e-01, 2.71358291e-06],
                           [9.95462789e-01, 2.80542716e-06],
                           [9.97632039e-01, 2.89832141e-06],
                           [9.93685586e-01, 2.99110166e-06],
                           [8.23528508e-01, 3.08399842e-06]] )

    assert_allclose( band, true_band, atol=1E-7 )


def setup_SF():
    ''' Initialize values for these parameters needed for the analyticnrm2 tests.

        Returns
        -------
        pixel (optional, via **kwargs) : float
            pixel scale

        fov : integer
            number of detector pixels on a side

        oversample : integer
            oversampling factor

        ctrs : float, float
            coordinates of hole centers

        d : float
            hole diameter

        lam : float
            wavelength

        phi : float
            distance of fringe from hole center in units of waves

        centering : string
            if set to 'PIXELCENTERED' or unspecified, the offsets will be set to
            (0.5,0.5); if set to 'PIXELCORNER', the offsets will be set to
            (0.0,0.0).
    '''
    pixel = 3.1125038327232215e-07
    fov = 2
    oversample = 3
    ctrs = np.array( [[ 0.06864653, -2.63910736],
                      [-2.28553695, -0.05944972],
                      [ 2.31986022, -1.26010406],
                      [-2.31986022,  1.26010406],
                      [-1.19424838,  1.94960579],
                      [ 2.25121368,  1.3790035 ],
                      [ 1.09127858,  2.00905525]]
                   )
    d = 0.8
    lam =  2.3965000082171173e-06
    phi = np.zeros(7, dtype=np.float32)
    centering = (0.5, 0.5)

    return pixel, fov, oversample, ctrs, d, lam, phi, centering


def setup_hexee():
    ''' Initialize values for parameters needed for the hexee tests.

        Returns
        -------
        xi : 2D float array
            hexagon's coordinate center at center of symmetry, along flat edge

        eta : 2D float array
            hexagon's coordinate center at center of symmetry, normal to xi;
            (not currently used)

        c (optional, via **kwargs) : tuple(float, float)
            coordinates of center

        pixel (optional, via **kwargs) : float
            pixel scale

        d (optional, via **kwargs) : float
            flat-to-flat distance across hexagon

        lam : (optional, via **kwargs) : float
            wavelength

        minus : (optional, via **kwargs) boolean
            if set, use flipped sign of xi in calculation
    '''
    xdim, ydim = 3, 3
    xi = np.zeros( ydim*xdim ).reshape(( ydim, xdim ))
    eta = np.zeros( ydim*xdim ).reshape(( ydim, xdim ))

    for ii in range(ydim):
        xi[ ii, : ] = ii
        eta[ :, ii] = ii

    kwargs =  {'d': 0.8, 'c': (28.0, 28.0), 'lam': 2.3965000082171173e-06,
                'pixel': 1.0375012775744072e-07, 'minus': False}

    return xi, eta, kwargs
