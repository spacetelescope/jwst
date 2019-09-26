"""

Unit tests for ami_analyze

"""
import numpy as np
import math

from jwst.ami import utils, leastsqnrm
from jwst.ami.leastsqnrm import hexpb, ffc, ffs, return_CAs
from jwst.ami.leastsqnrm import closurephase, redundant_cps
from jwst.ami.leastsqnrm import populate_symmamparray, populate_antisymmphasearray
from jwst.ami.leastsqnrm import tan2visibilities, model_array

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
    assert_allclose( hexpb_arr, true_hexpb_arr, atol=1E-8)

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
    assert_allclose( CAs, true_CAs, atol=1E-8 )


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

