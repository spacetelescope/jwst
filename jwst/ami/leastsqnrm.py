import logging

import numpy as np
import numpy.linalg as linalg
from scipy.special import comb, jv

from . import hexee


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def flip(holearray):
    """
    Short Summary
    -------------
    Change sign of 2nd coordinate of holes

    Parameters
    ----------
    holearray: 2D float array
        coordinates of holes

    Return
    ------
    fliparray: 2D float array
        flipped coordinates of holes
    """

    fliparray = holearray.copy()
    fliparray[:, 1] = -1 * holearray[:, 1]

    return fliparray


def rotatevectors(vectors, thetarad):
    """
    Short Summary
    -------------
    Rotate vectors by specified angle

    Parameters
    ----------
    vectors: 2D float array
        list of vectors - e.g. nrm hole centers; positive x decreases under
        slight rotation, and positive y increases under slight rotation

    thetarad: float
        rotation angle

    Returns
    -------
    rot_vectors: 2D float array
        rotated vectors
    """
    c, s = (np.cos(thetarad), np.sin(thetarad))
    ctrs_rotated = []
    for vector in vectors:
        ctrs_rotated.append([c * vector[0] - s * vector[1],
                             s * vector[0] + c * vector[1]])

    rot_vectors = np.array(ctrs_rotated)

    return rot_vectors


def mas2rad(mas):
    """
    Short Summary
    -------------
    Convert angle in milli arc-sec to radians

    Parameters
    ----------
    mas: float
        angle in milli arc-sec

    Returns
    -------
    rad: float
        angle in radians
    """

    rad = mas * (10**(-3)) / (3600 * 180 / np.pi)
    return rad


def rad2mas(rad):
    """
    Short Summary
    -------------
    Convert input angle in radians to milli arc sec

    Parameters
    ----------
    rad: float
        input angle in radians

    Returns
    -------
    mas: float
        input angle in milli arc sec
    """
    mas = rad * (3600. * 180 / np.pi) * 10.**3

    return mas


def sin2deltapistons(coeffs):
    """
    Short Summary
    -------------
    Each baseline has one sine and one cosine fringe with a coefficient that
    depends on the piston difference between the two holes that make the
    baseline.  For a 7-hole mask there are 21 baselines and therefore there
    are 42 sine and cosine terms that contribute to the fringe model. This
    function calculate the sine of this piston difference.

    Parameters
    ----------
    coeffs: 1D float array
        array of piston differences

    Returns
    -------
    delta: 1D float array
        sine of piston differences
    """
    asize = int((len(coeffs) -1)/2)

    delta = np.zeros(asize)
    for q in range(asize):
        delta[q] = np.arcsin(coeffs[2*q+2]) / (np.pi*2.0)

    return delta


def cos2deltapistons(coeffs):
    """
    Short Summary
    -------------
    Each baseline has one sine and one cosine fringe with a coefficient that
    depends on the piston difference between the two holes that make the
    baseline.  For a 7-hole mask there are 21 baselines and therefore there
    are 42 sine and cosine terms that contribute to the fringe model. This
    function calculate the cosine of this piston difference.

    Parameters
    ----------
    coeffs: 1D float array
        array of piston differences

    Returns
    -------
    delta: 1D float array
        cosine of piston differences
    """
    asize = int((len(coeffs) -1)/2)

    delta = np.zeros(asize)
    for q in range(asize):
        if coeffs[2*q+2]<0:
            sgn = -1
        else:
            sgn = 1
        delta[q] = sgn * np.arccos(coeffs[2*q+1]) / (np.pi*2.0)

    return delta


def replacenan(array):
    """
    Short Summary
    -------------
    Replace singularities encountered in the analytical hexagon Fourier
    transform with the analytically derived limits.

    Parameters
    ----------
    array: 2D float array
        input array

    Returns
    -------
    array: 2D float array
        input array with NaNs replaced with analytically derived limits
    """
    nanpos = np.where(np.isnan(array))
    array[nanpos] = np.pi / 4

    return array


def primarybeam(kx, ky):
    """
    Short Summary
    -------------
    Calculate the envelope intensity for circular holes & monochromatic light

    Parameters
    ----------
    kx, ky: float, float
        x-component and y-component of image plane (spatial frequency) vector

    Return
    ------
    env_int: 2D float array
        envelope intensity for circular holes & monochromatic light
    """
    R = (primarybeam.d / primarybeam.lam) * primarybeam.pitch *  \
            np.sqrt((kx - primarybeam.offx) * (kx - primarybeam.offx) + \
            (ky - primarybeam.offy) * (ky - primarybeam.offy))
    pb = replacenan(jv(1, np.pi * R) / (2.0 * R))

    pb = pb.transpose()

    env_int = pb * pb.conj()

    return env_int


def hexpb():
    """
    Short Summary
    -------------
    Calculate the primary beam for hexagonal holes.

    Parameters
    ----------
    None

    Returns
    -------
    pb * pb.conj(): 2D float array
        primary beam for hexagonal holes
    """
    pb = hexee.hex_eeAG(s=hexpb.size, c=(hexpb.offx, hexpb.offy),
            d=hexpb.d, lam=hexpb.lam, pitch=hexpb.pitch)

    return pb * pb.conj()


def ffc(kx, ky):
    """
    Short Summary
    -------------
    Calculate cosine terms of analytic model.

    Parameters
    ----------
    kx, ky: float, float
        x-component and y-component of image plane (spatial frequency) vector

    Returns
    -------
    cos_array: 2D float array
        cosine terms of analytic model
    """
    cos_array = 2 * np.cos(2 * np.pi * ffc.pitch *
                   ((kx - ffc.offx) * (ffc.ri[0] - ffc.rj[0]) +
                   (ky - ffc.offy) * (ffc.ri[1] - ffc.rj[1])) / ffc.lam)
    return cos_array


def ffs(kx, ky):
    """
    Short Summary
    -------------
    Calculate sine terms of analytic model.

    Parameters
    ----------
    kx, ky: float, float
        x-component and y-component of image plane (spatial frequency) vector

    Returns
    -------
    sin_array: 2D float array
        sine terms of analytic model
    """
    sin_array = -2 * np.sin(2 * np.pi * ffs.pitch *
                     ((kx - ffs.offx) * (ffs.ri[0] - ffs.rj[0]) +
                     (ky - ffs.offy) * (ffs.ri[1] - ffs.rj[1])) / ffs.lam)

    return sin_array


def model_array(ctrs, lam, oversample, pitch, fov, d,
                centering='PIXELCENTERED', shape='circ'):
    """
    Short Summary
    -------------
    Create a model using the specified wavelength.

    Parameters
    ----------
    ctrs: 2D float array
        centers of holes

    lam: float
        wavelength in the bandpass for this particular model

    oversample: integer
        oversampling factor

    pitch: float
        sampling pitch in radians in image plane

    fov: integer
        number of detector pixels on a side.

    d: float
        hole diameter for 'circ'; flat to flat distance for 'hex

    centering: string
        subpixel centering; for now only option is PIXELCENTERED, which means
        putting the brightest detector pixel at the center of the trimmed data
        frame or simulated image.

    shape: string
        shape of hole; possible values are 'circ', 'hex', and 'fringe'

    Returns
    -------
    if 'shape' == 'circ', returns the primary beam (2D float array)
        for circular holes.
    if 'shape' == 'hex', returns the primary beam (2D float array)
        for hexagonal holes.

    ffmodel: list of 3 2D float arrays
        model array
    """
    if centering == 'PIXELCORNER':
        off = np.array([0.0, 0.0])
    elif centering == 'PIXELCENTERED':
        off = np.array([0.5, 0.5])
    else:
        off = centering

    log.debug('------------------')
    log.debug('Model Parameters:')
    log.debug('------------------')
    log.debug('pitch:%s fov:%s oversampling:%s ', pitch, fov, oversample)
    log.debug('centers:%s', ctrs)
    log.debug('wavelength:%s  centering:%s off:%s ', lam, centering, off)
    log.debug('shape:%s d:%s ', shape, d)

    # primary beam parameters:
    primarybeam.shape = shape
    primarybeam.lam = lam
    primarybeam.size = (oversample * fov, oversample * fov)
    primarybeam.offx = oversample * fov / 2.0 - off[0] # in pixels
    primarybeam.offy = oversample * fov / 2.0 - off[1]
    primarybeam.pitch = pitch / float(oversample)
    primarybeam.d = d

    hexpb.shape = shape
    hexpb.lam = lam
    hexpb.size = (oversample * fov, oversample * fov)
    hexpb.offx = oversample * fov / 2.0 - off[0] # in pixels
    hexpb.offy = oversample * fov / 2.0 - off[1]
    hexpb.pitch = pitch / float(oversample)
    hexpb.d = d

    # model fringe matrix parameters:
    ffc.N = len(ctrs) # number of holes
    ffc.lam = lam
    ffc.over = oversample
    ffc.pitch = pitch / float(oversample)
    ffc.size = (oversample * fov, oversample * fov)
    ffc.offx = oversample * fov / 2.0 - off[0]
    ffc.offy = oversample * fov / 2.0 - off[1]

    ffs.N = len(ctrs) # number of holes
    ffs.lam = lam
    ffs.over = oversample
    ffs.pitch = pitch / float(oversample)
    ffs.size = (oversample * fov, oversample * fov)
    ffs.offx = oversample * fov / 2.0 - off[0]
    ffs.offy = oversample * fov / 2.0 - off[1]

    alist = []
    for i in range(ffc.N - 1):
        for j in range(ffc.N - 1):
            if j + i + 1 < ffc.N:
                alist = np.append(alist, i)
                alist = np.append(alist, j + i + 1)
    alist = alist.reshape(len(alist) // 2, 2)

    ffmodel = []
    ffmodel.append(ffc.N * np.ones(ffc.size))
    for q, r in enumerate(alist):
        # r[0] and r[1] are holes i and j, x-coord: 0, y-coord: 1
        ffc.ri = ctrs[int(r[0])]
        ffc.rj = ctrs[int(r[1])]
        ffs.ri = ctrs[int(r[0])]
        ffs.rj = ctrs[int(r[1])]

        ffmodel.append( np.transpose(np.fromfunction(ffc, ffc.size)) )
        ffmodel.append( np.transpose(np.fromfunction(ffs, ffs.size)) )

    if shape == 'circ': # if unspecified (default), or specified as 'circ'
        return np.fromfunction(primarybeam, ffc.size), ffmodel
    elif shape == 'hex':
        return hexpb(), ffmodel
    else:
        log.critical('Must provide a valid hole shape. Current supported shapes \
        are circ and hex.')
        return None


def weighted_operations(img, model, weights):
    """
    Short Summary
    -------------
    Performs least squares matrix operations to solve A x = b, where A is the
    model, b is the data (image), and x is the coefficient vector we are solving
    for.  In 2-D, data x = inv(At.A).(At.b)

    Parameters
    ----------
    img: 2D float array
        input data

    model: 2D float array
        analytic model

    weights: 2D float array
        values of weights

    Returns
    -------
    x: 1D float array
        coefficient vector

    res: 2D float array
        residual; difference between model and fit
    """
    clist = weights.reshape(weights.shape[0] * weights.shape[1])**2
    flatimg = img.reshape(np.shape(img)[0] * np.shape(img)[1])
    nanlist = np.where(np.isnan(flatimg))
    flatimg = np.delete(flatimg, nanlist)
    clist = np.delete(clist, nanlist)
    # A
    flatmodel_nan = model.reshape(np.shape(model)[0] * np.shape(model)[1],
                    np.shape(model)[2])
    flatmodel = np.zeros((len(flatimg), np.shape(model)[2]))

    for fringe in range(np.shape(model)[2]):
        flatmodel[:,fringe] = np.delete(flatmodel_nan[:,fringe], nanlist)

    # At (A transpose)
    flatmodeltransp = flatmodel.transpose()
    # At.C.A (makes square matrix)
    CdotA = flatmodel.copy()

    for i in range(flatmodel.shape[1]):
        CdotA[:,i] = clist * flatmodel[:,i]

    modelproduct = np.dot(flatmodeltransp, CdotA)
    # At.C.b
    Cdotb = clist * flatimg
    data_vector = np.dot(flatmodeltransp, Cdotb)
    # inv(At.C.A)
    inverse = linalg.inv(modelproduct)

    x = np.dot(inverse, data_vector)
    res = np.dot(flatmodel, x) - flatimg
    naninsert = nanlist[0] - np.arange(len(nanlist[0]))
    res = np.insert(res, naninsert, np.nan)
    res = res.reshape(img.shape[0], img.shape[1])

    return x, res


def matrix_operations(img, model, flux=None, linfit=False):
    """
    Short Summary
    -------------
    Use least squares matrix operations to solve A x = b, where A is the model,
    b is the data (img), and x is the coefficient vector we are solving for.
    In 2-D, data x = inv(At.A).(At.b).  If a flux is given, it will be used it
    to normalize the data.

    Parameters
    ----------
    img: 2D float array
        input data

    model: 2D float array
        analytic model

    flux: float
        normalization factor

    Returns
    -------
    x: 1D float array
        solution to fit

    res: 2D float array
        residuals in fit

    cond: float
        condition number of the inverse of the product of model and its
        transpose
    """
    flatimg = img.reshape(np.shape(img)[0] * np.shape(img)[1])
    nanlist = np.where(np.isnan(flatimg))
    flatimg = np.delete(flatimg, nanlist)

    if flux is not None:
        flatimg = flux * flatimg / flatimg.sum()

    # A
    flatmodel_nan = model.reshape(np.shape(model)[0] * np.shape(model)[1],
                          np.shape(model)[2])

    flatmodel = np.zeros((len(flatimg), np.shape(model)[2]))

    log.debug('Matrix_opers - flat model dimensions: %s', np.shape(flatmodel))
    log.debug('Matrix_opers - flat image dimensions: %s', np.shape(flatimg))

    for fringe in range(np.shape(model)[2]):
        flatmodel[:, fringe] = np.delete(flatmodel_nan[:, fringe], nanlist)

    # At (A transpose)
    flatmodeltransp = flatmodel.transpose()
    # At.A (makes square matrix)
    modelproduct = np.dot(flatmodeltransp, flatmodel)
    # At.b
    data_vector = np.dot(flatmodeltransp, flatimg)
    # inv(At.A)
    inverse = linalg.inv(modelproduct)
    cond = np.linalg.cond(inverse)

    x = np.dot(inverse, data_vector)
    res = np.dot(flatmodel, x) - flatimg
    naninsert = nanlist[0] - np.arange(len(nanlist[0]))
    res = np.insert(res, naninsert, np.nan)
    res = res.reshape(img.shape[0], img.shape[1])

    log.debug('------------------')
    log.debug('Matrix Operations:')
    log.debug('------------------')
    log.debug('model flux:%s data flux:%s flat model dimensions:%s ', flux,
              flatimg.sum(), np.shape(flatmodel))
    log.debug('flat image dimensions:%s model transpose dimensions:%s ',
              np.shape(flatimg), np.shape(flatmodeltransp))
    log.debug('transpose * image data dimensions:%s flatimg * transpose' +
              'dimensions:%s ', np.shape(data_vector), np.shape(inverse))

    if linfit:
        try:
            from linearfit import linearfit
            # dependent variables
            M = np.mat(flatimg)

            # photon noise
            noise = np.sqrt(np.abs(flatimg))

            # this sets the weights of pixels fulfilling condition to zero
            weights = np.where(np.abs(flatimg)<=1.0, 0.0, 1.0/(noise**2))

            # uniform weight
            wy = weights
            S = np.mat(np.diag(wy));
            # matrix of independent variables
            C = np.mat(flatmodeltransp)

            # initialize object
            result = linearfit.LinearFit(M,S,C)

            # do the fit
            result.fit()

            # delete inverse_covariance_matrix to reduce size of pickled file
            result.inverse_covariance_matrix = []

            linfit_result = result

        except ImportError:
            linfit_result = None
    else:
        linfit_result = None

    return x, res, cond, linfit_result


def multiplyenv(env, fringeterms):
    """
    Short Summary
    -------------
    Multiply the envelope by each fringe 'image'.

    Parameters
    ----------
    env: 2D float array
        envelope

    fringeterms: list of 3 2D float arrays
        model

    Returns
    -------
    full: 3D float array
        envelope multiplied by each fringe 'image'
    """
    # The envelope has size (fov, fov). This multiplies the envelope by each
    #    of the 43 slices in the fringe model
    full = np.ones((np.shape(fringeterms)[1], np.shape(fringeterms)[2],
                    np.shape(fringeterms)[0] + 1))

    for i, val in enumerate(fringeterms):
        full[:, :, i] = env * fringeterms[i]

    log.debug('Total number of fringe terms: %s', len(fringeterms) - 1)

    return full


def tan2visibilities(coeffs):
    """
    Long Summary
    ------------
    Technically the fit measures phase AND amplitude, so to retrieve the
    phase we need to consider both sin and cos terms. Consider one fringe:
    A { cos(kx)cos(dphi) + sin(kx)sin(dphi) } =
    A(a cos(kx) + b sin(kx)), where a = cos(dphi) and b = sin(dphi)
    and A is the fringe amplitude, therefore coupling a and b.
    In practice we measure A*a and A*b from the coefficients, so:
    Ab/Aa = b/a = tan(dphi)
    call a' = A*a and b' = A*b (we actually measure a', b')
    (A*sin(dphi))^2 + (A*cos(dphi)^2) = A^2 = a'^2 + b'^2

    Short Summary
    -------------
    From the solution to the fit, calculate the fringe amplitude and phase.

    Parameters
    ----------
    coeffs: 1D float array

    Returns
    -------
    amp, delta: 1D float array, 1D float array
        fringe amplitude & phase
    """
    delta = np.zeros(int((len(coeffs) -1)/2))
    amp = np.zeros(int((len(coeffs) -1)/2))
    for q in range(int((len(coeffs) -1)/2)):
        delta[q] = (np.arctan2(coeffs[2*q+2], coeffs[2*q+1]))
        amp[q] = np.sqrt(coeffs[2*q+2]**2 + coeffs[2*q+1]**2)

    log.debug(f"tan2visibilities: shape coeffs:{np.shape(coeffs)} "
        f"shape delta:{np.shape(delta)}")

    # returns fringe amplitude & phase
    return amp, delta


def populate_antisymmphasearray(deltaps, n=7):
    """
    Short Summary
    -------------
    Populate the antisymmetric fringe phase array:

    fringephasearray[0,q+1:] = coeffs[0:6]
    fringephasearray[1,q+2:] = coeffs[6:11]
    fringephasearray[2,q+3:] = coeffs[11:15]
    fringephasearray[3,q+4:] = coeffs[15:18]
    fringephasearray[4,q+5:] = coeffs[18:20]
    fringephasearray[5,q+6:] = coeffs[20:]

    Parameters
    ----------
    deltaps: 1D float array
        pistons between each pair of holes

    n: integer
        number of holes

    Returns
    -------
    arr: 2D float array
        fringe phases between each pair of holes
    """
    # Initialize fringe phase array
    arr = np.zeros((n, n))

    step = 0
    n = n - 1
    for h in range(n):
        arr[h, h+1:] = deltaps[step:step+n]
        step += n
        n -= 1

    arr -= arr.T

    return arr


def populate_symmamparray(amps, n=7):
    """
    Short Summary
    -------------
    Populate the symmetric fringe amplitude array

    Parameters
    ----------
    amps: 1D float array
        fringe visibility between each pair of holes

    n: integer
        number of holes

    Returns
    -------
    arr: 2D float array
        fringe amplitude array
    """
    arr = np.zeros((n, n))

    step = 0
    n = n - 1

    for h in range(n):
        arr[h, h + 1:] = amps[step:step + n]
        step += n
        n -= 1

    arr += arr.T

    return arr


def redundant_cps(deltaps, n=7):
    """
    Short Summary
    -------------
    Calculate closure phases for each set of 3 holes

    Parameters
    ----------
    deltaps: 1D float array
        pistons between each pair of holes

    n: integer
        number of holes

    Returns
    -------
    cps: 1D float array
        closure phases
    """
    arr = populate_antisymmphasearray(deltaps, n=n) # fringe phase array

    cps = np.zeros(int(comb(n, 3)))

    nn = 0
    for kk in range(n - 2):
        for ii in range(n - kk - 2):
            for jj in range(n - kk - ii - 2):
                cps[nn + jj] = arr[kk, ii + kk + 1] \
                       + arr[ii + kk + 1, jj + ii + kk + 2] \
                       + arr[jj + ii + kk + 2, kk]

            nn += jj + 1

    return cps


def closurephase(deltap, n=7):
    """
    Short Summary
    -------------
    Calculate closure phases between each pair of holes

    Parameters
    ----------
    deltap: 1D float array
        pistons between each pair of holes

    n: integer
        number of holes in the mask; 7 and 10 holes available (JWST & GPI))

    Returns
    -------
    cps: 1D float array
        closure phases
    """
    # p is a triangular matrix set up to calculate closure phases
    if n == 7:
        p = np.array([deltap[:6], deltap[6:11], deltap[11:15],
                deltap[15:18], deltap[18:20], deltap[20:]], dtype=object)
    elif n == 10:
        p = np.array([deltap[:9], deltap[9:17], deltap[17:24],
                deltap[24:30], deltap[30:35], deltap[35:39],
                deltap[39:42], deltap[42:44], deltap[44:]], dtype=object)
    else:
        log.critical('invalid hole number: %s', n)

    # calculates closure phases for general N-hole mask (with p-array set
    #     up properly above)
    cps = np.zeros((n - 1) * (n - 2) // 2)
    for j1 in range(n - 2):
        for j2 in range(n - 2 - j1):
            cps[int(j1 * ((n + (n - 3) - j1) / 2.0)) + j2] = \
                p[j1][0] + p[j1 + 1][j2] - p[j1][j2 + 1]

    return cps


def closure_amplitudes(amps, n=7):
    """
    Short Summary
    -------------
    Calculate closure amplitudes

    Parameters
    ----------
    amps: 1D float array
         fringe amplitudes

    N: integer
        number of holes

    Returns
    -------
    CAs: 1D float array
        closure amplitudes
    """
    arr = populate_symmamparray(amps, n=n) # fringe amp array
    nn = 0

    cas = np.zeros(int(comb(n, 4)))

    for ii in range(n - 3):
        for jj in range(n - ii - 3):
            for kk in range(n - jj - ii - 3):
                for ll in range(n - jj - ii - kk - 3):
                    cas[nn + ll] = arr[ii, jj + ii + 1] \
                           * arr[ll + ii + jj + kk + 3, kk + jj + ii + 2] \
                           / (arr[ii, kk + ii + jj + 2] * \
                              arr[jj + ii + 1, ll + ii + jj + kk + 3])
                nn = nn + ll + 1

    return cas
