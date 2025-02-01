import logging

import numpy as np
import numpy.linalg as linalg
from scipy.special import comb, jv

from . import hexee


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def flip(holearray):
    """
    Change sign of 2nd coordinate of holes.

    Parameters
    ----------
    holearray : 2D float array
        Coordinates of holes

    Returns
    -------
    fliparray : 2D float array
        Flipped coordinates of holes
    """
    fliparray = holearray.copy()
    fliparray[:, 1] = -1 * holearray[:, 1]

    return fliparray


def rotatevectors(vectors, thetarad):
    """
    Rotate vectors by specified angle.

    Parameters
    ----------
    vectors : 2D float array
        List of vectors - e.g. nrm hole centers; positive x decreases under
        slight rotation, and positive y increases under slight rotation

    thetarad : float
        Rotation angle

    Returns
    -------
    rot_vectors : 2D float array
        Rotated vectors
    """
    c, s = (np.cos(thetarad), np.sin(thetarad))
    ctrs_rotated = []
    for vector in vectors:
        ctrs_rotated.append([c * vector[0] - s * vector[1], s * vector[0] + c * vector[1]])

    rot_vectors = np.array(ctrs_rotated)

    return rot_vectors


def mas2rad(mas):
    """
    Convert angle in milli arc-sec to radians.

    Parameters
    ----------
    mas : float
        Angle in milli arc-sec

    Returns
    -------
    rad : float
        Angle in radians
    """
    rad = mas * (10 ** (-3)) / (3600 * 180 / np.pi)
    return rad


def rad2mas(rad):
    """
    Convert input angle in radians to milli arc sec.

    Parameters
    ----------
    rad : float
        Input angle in radians

    Returns
    -------
    mas : float
        Input angle in milli arc sec
    """
    mas = rad * (3600.0 * 180 / np.pi) * 10.0**3

    return mas


def sin2deltapistons(coeffs):
    """
    Calculate the sine of the piston difference.

    Each baseline has one sine and one cosine fringe with a coefficient that
    depends on the piston difference between the two holes that make the
    baseline.  For a 7-hole mask there are 21 baselines and therefore there
    are 42 sine and cosine terms that contribute to the fringe model. This
    function calculate the sine of this piston difference.

    Parameters
    ----------
    coeffs : 1D float array
        Array of piston differences

    Returns
    -------
    delta : 1D float array
        Sine of piston differences
    """
    asize = int((len(coeffs) - 1) / 2)

    delta = np.zeros(asize)
    for q in range(asize):
        delta[q] = np.arcsin(coeffs[2 * q + 2]) / (np.pi * 2.0)

    return delta


def cos2deltapistons(coeffs):
    """
    Calculate the cosine of the piston difference.

    Each baseline has one sine and one cosine fringe with a coefficient that
    depends on the piston difference between the two holes that make the
    baseline.  For a 7-hole mask there are 21 baselines and therefore there
    are 42 sine and cosine terms that contribute to the fringe model. This
    function calculates the cosine of this piston difference.

    Parameters
    ----------
    coeffs : 1D float array
        Array of piston differences

    Returns
    -------
    delta : 1D float array
        Cosine of piston differences
    """
    asize = int((len(coeffs) - 1) / 2)

    delta = np.zeros(asize)
    for q in range(asize):
        if coeffs[2 * q + 2] < 0:
            sgn = -1
        else:
            sgn = 1
        delta[q] = sgn * np.arccos(coeffs[2 * q + 1]) / (np.pi * 2.0)

    return delta


def replacenan(array):
    """
    Replace singularities in analytical hexagon Fourier transform.

    Replace NaN values with the analytically derived limits.

    Parameters
    ----------
    array : 2D float array
        Input array

    Returns
    -------
    array : 2D float array
        Input array with NaNs replaced with analytically derived limits
    """
    nanpos = np.where(np.isnan(array))
    array[nanpos] = np.pi / 4

    return array


def primarybeam(kx, ky):
    """
    Calculate the envelope intensity for circular holes & monochromatic light.

    Parameters
    ----------
    kx, ky : float, float
        The x-component and y-component of image plane (spatial frequency) vector

    Returns
    -------
    env_int : 2D float array
        Envelope intensity for circular holes & monochromatic light
    """
    r = (
        (primarybeam.d / primarybeam.lam)
        * primarybeam.pitch
        * np.sqrt(
            (kx - primarybeam.offx) * (kx - primarybeam.offx)
            + (ky - primarybeam.offy) * (ky - primarybeam.offy)
        )
    )
    pb = replacenan(jv(1, np.pi * r) / (2.0 * r))

    pb = pb.transpose()

    env_int = pb * pb.conj()

    return env_int


def hexpb():
    """
    Calculate the primary beam for hexagonal holes.

    Returns
    -------
    pb * pb.conj() : 2D float array
        Primary beam for hexagonal holes
    """
    pb = hexee.hex_eeag(
        s=hexpb.size,
        c=(hexpb.offx, hexpb.offy),
        d=hexpb.d,
        lam=hexpb.lam,
        pitch=hexpb.pitch,
    )

    return pb * pb.conj()


def ffc(kx, ky):
    """
    Calculate cosine terms of analytic model.

    Parameters
    ----------
    kx, ky : float, float
        The x-component and y-component of image plane (spatial frequency) vector

    Returns
    -------
    cos_array : 2D float array
        Cosine terms of analytic model
    """
    cos_array = 2 * np.cos(
        2
        * np.pi
        * ffc.pitch
        * ((kx - ffc.offx) * (ffc.ri[0] - ffc.rj[0]) + (ky - ffc.offy) * (ffc.ri[1] - ffc.rj[1]))
        / ffc.lam
    )
    return cos_array


def ffs(kx, ky):
    """
    Calculate sine terms of analytic model.

    Parameters
    ----------
    kx, ky : float, float
        The x-component and y-component of image plane (spatial frequency) vector

    Returns
    -------
    sin_array : 2D float array
        Sine terms of analytic model
    """
    sin_array = -2 * np.sin(
        2
        * np.pi
        * ffs.pitch
        * ((kx - ffs.offx) * (ffs.ri[0] - ffs.rj[0]) + (ky - ffs.offy) * (ffs.ri[1] - ffs.rj[1]))
        / ffs.lam
    )

    return sin_array


def model_array(ctrs, lam, oversample, pitch, fov, d, centering="PIXELCENTERED", shape="circ"):
    """
    Create a model using the specified wavelength.

    Parameters
    ----------
    ctrs : 2D float array
        Centers of holes

    lam : float
        Wavelength in the bandpass for this particular model

    oversample : int
        Oversampling factor

    pitch : float
        Sampling pitch in radians in image plane

    fov : int
        Number of detector pixels on a side.

    d : float
        Hole diameter for 'circ'; flat to flat distance for 'hex

    centering : str
        Subpixel centering; for now only option is PIXELCENTERED, which means
        putting the brightest detector pixel at the center of the trimmed data
        frame or simulated image.

    shape : str
        Shape of hole; possible values are 'circ', 'hex', and 'fringe'

    Returns
    -------
    ffmodel : list of 3 2D float arrays
        Model array
        if 'shape' == 'circ', returns the primary beam (2D float array)
            for circular holes.
        if 'shape' == 'hex', returns the primary beam (2D float array)
            for hexagonal holes.
    """
    if centering == "PIXELCORNER":
        off = np.array([0.0, 0.0])
    elif centering == "PIXELCENTERED":
        off = np.array([0.5, 0.5])
    else:
        off = centering

    log.debug("------------------")
    log.debug("Model Parameters:")
    log.debug("------------------")
    log.debug("pitch:%s fov:%s oversampling:%s ", pitch, fov, oversample)
    log.debug("centers:%s", ctrs)
    log.debug("wavelength:%s  centering:%s off:%s ", lam, centering, off)
    log.debug("shape:%s d:%s ", shape, d)

    # primary beam parameters:
    primarybeam.shape = shape
    primarybeam.lam = lam
    primarybeam.size = (oversample * fov, oversample * fov)
    primarybeam.offx = oversample * fov / 2.0 - off[0]  # in pixels
    primarybeam.offy = oversample * fov / 2.0 - off[1]
    primarybeam.pitch = pitch / float(oversample)
    primarybeam.d = d

    hexpb.shape = shape
    hexpb.lam = lam
    hexpb.size = (oversample * fov, oversample * fov)
    hexpb.offx = oversample * fov / 2.0 - off[0]  # in pixels
    hexpb.offy = oversample * fov / 2.0 - off[1]
    hexpb.pitch = pitch / float(oversample)
    hexpb.d = d

    # model fringe matrix parameters:
    ffc.N = len(ctrs)  # number of holes
    ffc.lam = lam
    ffc.over = oversample
    ffc.pitch = pitch / float(oversample)
    ffc.size = (oversample * fov, oversample * fov)
    ffc.offx = oversample * fov / 2.0 - off[0]
    ffc.offy = oversample * fov / 2.0 - off[1]

    ffs.N = len(ctrs)  # number of holes
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
    for r in alist:
        # r[0] and r[1] are holes i and j, x-coord: 0, y-coord: 1
        ffc.ri = ctrs[int(r[0])]
        ffc.rj = ctrs[int(r[1])]
        ffs.ri = ctrs[int(r[0])]
        ffs.rj = ctrs[int(r[1])]

        ffmodel.append(np.transpose(np.fromfunction(ffc, ffc.size)))
        ffmodel.append(np.transpose(np.fromfunction(ffs, ffs.size)))

    if shape == "circ":  # if unspecified (default), or specified as 'circ'
        return np.fromfunction(primarybeam, ffc.size), ffmodel
    elif shape == "hex":
        return hexpb(), ffmodel
    else:
        log.critical(
            "Must provide a valid hole shape. Current supported shapes \
        are circ and hex."
        )
        return None


def weighted_operations(img, model, dqm=None):
    """
    Perform least squares matrix operations to solve A x = b weighting by Poisson variance.

    A is the model, b is the data (image), and x is the coefficient vector we are solving for.

    Here we are weighting data by Poisson variance:
      x = inv(At.W.A).(At.W.b)
      where W is a diagonal matrix of weights w_i,
      weighting each data point i by the inverse of its variance:
         w_i = 1 / sigma_i^2
      For photon noise, the data, i.e. the image values b_i  have variance
      proportional to b_i with an e.g. ADU to electrons conversion factor.
      If this factor is the same for all pixels, we do not need to include
      it here.

    Parameters
    ----------
    img : 2D float array
        Input data

    model : 2D float array
        Analytic model

    dqm : 2D bool array
        Bad pixel mask

    Returns
    -------
    x : 1D float array
        Coefficient vector

    res : 2D float array
        Residual; difference between model and fit

    Notes
    -----
    Use matrix_operations() for equal weighting of data.
    """
    # Remove not-to-be-fit data from the flattened "img" data vector
    flatimg = img.reshape(np.shape(img)[0] * np.shape(img)[1])
    flatdqm = dqm.reshape(np.shape(img)[0] * np.shape(img)[1])

    if dqm is not None:
        nanlist = np.where(flatdqm)  # where DO_NOT_USE up.
    else:
        nanlist = (np.array(()),)  # shouldn't occur w/MAST JWST data

    # see original linearfit https://github.com/agreenbaum/ImPlaneIA:
    # agreenbaum committed on May 21, 2017 1 parent 3e0fb8b
    # commit bf02eb52c5813cb5d77036174a7caba703f9d366
    #
    flatimg = np.delete(flatimg, nanlist)  # DATA values

    # photon noise variance - proportional to ADU
    # (for roughly uniform adu2electron factor)
    variance = np.abs(flatimg)
    # this resets the weights of pixels with negative or unity values to zero
    # we ignore data with unity or lower values - weight it not-at-all..
    weights = np.where(flatimg <= 1.0, 0.0, 1.0 / np.sqrt(variance))  # anand 2022 Jan

    log.debug(f"{len(nanlist[0]):d} bad pixels skipped in weighted fringefitter")

    # A - but delete all pixels flagged by dq array
    flatmodel_nan = model.reshape(np.shape(model)[0] * np.shape(model)[1], np.shape(model)[2])
    flatmodel = np.zeros((len(flatimg), np.shape(model)[2]))
    for fringe in range(np.shape(model)[2]):
        flatmodel[:, fringe] = np.delete(flatmodel_nan[:, fringe], nanlist)

    # A.w
    aw = flatmodel * weights[:, np.newaxis]
    bw = flatimg * weights
    # resids are pixel value residuals, flattened to 1d vector
    x, _rss, _rank, singvals = np.linalg.lstsq(aw, bw, rcond=None)

    # actual residuals in image:
    res = flatimg - np.dot(flatmodel, x)

    # put bad pixels back
    naninsert = nanlist[0] - np.arange(len(nanlist[0]))
    # calculate residuals with fixed but unused bad pixels as nans
    res = np.insert(res, naninsert, np.nan)
    res = res.reshape(img.shape[0], img.shape[1])

    cond = None
    return x, res, cond, singvals  # no condition number yet...


def matrix_operations(img, model, flux=None, linfit=False, dqm=None):
    """
    Use least squares matrix operations to solve A x = b.

    A is the model, b is the data (img), and x is the coefficient vector we are solving for.
    In 2-D, data x = inv(At.A).(At.b).  If a flux is given, it will be used it
    to normalize the data.

    Parameters
    ----------
    img : 2D float array
        Input data

    model : 2D float array
        Analytic model

    flux : float
        Normalization factor

    linfit : bool
        Whether to perform linear fit

    dqm : 2D bool array
        Bad pixel mask slice

    Returns
    -------
    x : 1D float array
        Solution to fit

    res : 2D float array
        Residuals in fit

    cond : float
        Condition number of the inverse of the product of model and its
        transpose
    """
    flatimg = img.reshape(np.shape(img)[0] * np.shape(img)[1])
    flatdqm = dqm.reshape(np.shape(img)[0] * np.shape(img)[1])
    log.info("fringefitting.leastsqnrm.matrix_operations(): ")
    log.info(f"\timg {img.shape:}")
    log.info(f"\tdqm {dqm.shape:}")
    log.info(
        f"\tL x W = {img.shape[0]:d} x {img.shape[1]:d} = {img.shape[0] * img.shape[1]:d}",
    )
    log.info(f"\tflatimg {flatimg.shape:}")
    log.info(f"\tflatdqm {flatdqm.shape:}")

    log.info("")
    log.info("\ttype(dqm) %s", type(dqm))
    if dqm is not None:
        nanlist = np.where(flatdqm)  # where DO_NOT_USE up.
    else:
        nanlist = (np.array(()),)  # shouldn't occur w/MAST JWST data

    log.info(f"\ttype(nanlist) {type(nanlist):}, len={len(nanlist):}")
    log.info(f"\tnumber of nanlist pixels: {len(nanlist[0]):d} items")
    log.info(f"\t{len(nanlist[0]):d} DO_NOT_USE pixels found in data slice")

    flatimg = np.delete(flatimg, nanlist)

    log.info(f"\tflatimg {flatimg.shape:} after deleting {len(nanlist[0]):d}")

    if flux is not None:
        flatimg = flux * flatimg / flatimg.sum()

    # A
    flatmodel_nan = model.reshape(np.shape(model)[0] * np.shape(model)[1], np.shape(model)[2])
    flatmodel = np.zeros((len(flatimg), np.shape(model)[2]))
    log.info(f"\tflatmodel_nan {flatmodel_nan.shape:}")
    log.info(f"\tflatmodel     {flatmodel.shape:}")
    log.info(f"\tdifference    {flatmodel_nan.shape[0] - flatmodel.shape[0]:}")
    log.info("flat model dimensions %s", np.shape(flatmodel))
    log.info("flat image dimensions %s", np.shape(flatimg))

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
    res = flatimg - np.dot(flatmodel, x)

    # put bad pixels back
    naninsert = nanlist[0] - np.arange(len(nanlist[0]))
    # calculate residuals with fixed but unused bad pixels as nans
    res = np.insert(res, naninsert, np.nan)
    res = res.reshape(img.shape[0], img.shape[1])

    log.info("model flux %s", flux)
    log.info("data flux %s", flatimg.sum())
    log.info("flat model dimensions %s", np.shape(flatmodel))
    log.info("model transpose dimensions %s", np.shape(flatmodeltransp))
    log.info("flat image dimensions %s", np.shape(flatimg))
    log.info("transpose * image data dimensions %s", np.shape(data_vector))
    log.info("flat img * transpose dimensions %s", np.shape(inverse))

    if linfit:
        try:
            from linearfit import linearfit  # type: ignore[import-not-found]

            # dependent variables
            M = np.asmatrix(flatimg)  # noqa: N806

            # photon noise
            noise = np.sqrt(np.abs(flatimg))

            # this sets the weights of pixels fulfilling condition to zero
            weights = np.where(np.abs(flatimg) <= 1.0, 0.0, 1.0 / (noise**2))

            # uniform weight
            wy = weights
            S = np.asmatrix(np.diag(wy))  # noqa: N806
            # matrix of independent variables
            C = np.asmatrix(flatmodeltransp)  # noqa: N806

            # initialize object
            result = linearfit.LinearFit(M, S, C)

            # do the fit
            result.fit()

            # delete inverse_covariance_matrix to reduce size of pickled file
            result.inverse_covariance_matrix = []

            linfit_result = result
            log.info("Returned linearfit result")

        except ImportError:
            linfit_result = None
            log.info("linearfit module not imported, no covariances saved.")
    else:
        linfit_result = None
        log.info("linearfit not attempted, no covariances saved.")

    return x, res, cond, linfit_result


def multiplyenv(env, fringeterms):
    """
    Multiply the envelope by each fringe 'image'.

    Parameters
    ----------
    env : 2D float array
        Envelope

    fringeterms : list of 3 2D float arrays
        Model

    Returns
    -------
    full : 3D float array
        Envelope multiplied by each fringe 'image'
    """
    # The envelope has size (fov, fov). This multiplies the envelope by each
    #    of the 43 slices in the fringe model
    full = np.ones(
        (
            np.shape(fringeterms)[1],
            np.shape(fringeterms)[2],
            np.shape(fringeterms)[0] + 1,
        )
    )

    for i, fringeterm in enumerate(fringeterms):
        full[:, :, i] = env * fringeterm

    log.debug("Total number of fringe terms: %s", len(fringeterms) - 1)

    return full


def tan2visibilities(coeffs):
    """
    From the solution to the fit, calculate the fringe amplitude and phase.

    Parameters
    ----------
    coeffs : 1D float array
        The coefficients providing the fit solution.

    Returns
    -------
    amp, delta : 1D float array, 1D float array
        Fringe amplitude & phase

    Notes
    -----
    Technically the fit measures phase AND amplitude, so to retrieve the
    phase we need to consider both sin and cos terms. Consider one fringe:
    A { cos(kx)cos(dphi) + sin(kx)sin(dphi) } =
    A(a cos(kx) + b sin(kx)), where a = cos(dphi) and b = sin(dphi)
    and A is the fringe amplitude, therefore coupling a and b.
    In practice we measure A*a and A*b from the coefficients, so:
    Ab/Aa = b/a = tan(dphi)
    call a' = A*a and b' = A*b (we actually measure a', b')
    (A*sin(dphi))^2 + (A*cos(dphi)^2) = A^2 = a'^2 + b'^2
    """
    delta = np.zeros(int((len(coeffs) - 1) / 2))
    amp = np.zeros(int((len(coeffs) - 1) / 2))
    for q in range(int((len(coeffs) - 1) / 2)):
        delta[q] = np.arctan2(coeffs[2 * q + 2], coeffs[2 * q + 1])
        amp[q] = np.sqrt(coeffs[2 * q + 2] ** 2 + coeffs[2 * q + 1] ** 2)

    log.debug(f"tan2visibilities: shape coeffs:{np.shape(coeffs)} shape delta:{np.shape(delta)}")

    # returns fringe amplitude & phase
    return amp, delta


def populate_antisymmphasearray(deltaps, n=7):
    """
    Populate the antisymmetric fringe phase array:.

    This array takes the form:

    fringephasearray[0,q+1:] = coeffs[0:6]
    fringephasearray[1,q+2:] = coeffs[6:11]
    fringephasearray[2,q+3:] = coeffs[11:15]
    fringephasearray[3,q+4:] = coeffs[15:18]
    fringephasearray[4,q+5:] = coeffs[18:20]
    fringephasearray[5,q+6:] = coeffs[20:]

    Parameters
    ----------
    deltaps : 1D float array
        Pistons between each pair of holes

    n : int, optional
        Number of holes (default=7)

    Returns
    -------
    arr : 2D float array
        Fringe phases between each pair of holes
    """
    # Initialize fringe phase array
    arr = np.zeros((n, n))

    step = 0
    n = n - 1
    for h in range(n):
        arr[h, h + 1 :] = deltaps[step : step + n]
        step += n
        n -= 1

    arr -= arr.T

    return arr


def populate_symmamparray(amps, n=7):
    """
    Populate the symmetric fringe amplitude array.

    Parameters
    ----------
    amps : 1D float array
        Fringe visibility between each pair of holes

    n : int, optional
        Number of holes (default=7)

    Returns
    -------
    arr: 2D float array
        Fringe amplitude array
    """
    arr = np.zeros((n, n))

    step = 0
    n = n - 1

    for h in range(n):
        arr[h, h + 1 :] = amps[step : step + n]
        step += n
        n -= 1

    arr += arr.T

    return arr


def t3_amplitudes(amps, n=7):
    """
    Populate the triple-product amplitude array (NOT closure amplitudes).

    Parameters
    ----------
    amps : 1D float array
        Fringe visibility between each pair of holes

    n : int, optional
        Number of holes (default=7)

    Returns
    -------
    cpamps : 1D float array
        Triple product amplitude array
    """
    arr = populate_symmamparray(amps, n=n)

    cpamps = np.zeros(int(comb(n, 3)))

    nn = 0
    for kk in range(n - 2):
        for ii in range(n - kk - 2):
            for jj in range(n - kk - ii - 2):
                cpamps[nn + jj] = (
                    arr[kk, ii + kk + 1]
                    * arr[ii + kk + 1, jj + ii + kk + 2]
                    * arr[jj + ii + kk + 2, kk]
                )

            nn += jj + 1

    return cpamps


def redundant_cps(deltaps, n=7):
    """
    Calculate closure phases for each set of 3 holes.

    Parameters
    ----------
    deltaps : 1D float array
        Pistons between each pair of holes

    n : int, optional
        Number of holes (default=7)

    Returns
    -------
    cps : 1D float array
        Closure phases
    """
    arr = populate_antisymmphasearray(deltaps, n=n)  # fringe phase array

    cps = np.zeros(int(comb(n, 3)))

    nn = 0
    for kk in range(n - 2):
        for ii in range(n - kk - 2):
            for jj in range(n - kk - ii - 2):
                cps[nn + jj] = (
                    arr[kk, ii + kk + 1]
                    + arr[ii + kk + 1, jj + ii + kk + 2]
                    + arr[jj + ii + kk + 2, kk]
                )

            nn += jj + 1

    return cps


def closurephase(deltap, n=7):
    """
    Calculate closure phases between each pair of holes.

    Parameters
    ----------
    deltap : 1D float array
        Pistons between each pair of holes

    n : int
        Number of holes in the mask; 7 and 10 holes available (JWST & GPI))

    Returns
    -------
    cps : 1D float array
        Closure phases
    """
    # p is a triangular matrix set up to calculate closure phases
    if n == 7:
        p = np.array(
            [
                deltap[:6],
                deltap[6:11],
                deltap[11:15],
                deltap[15:18],
                deltap[18:20],
                deltap[20:],
            ],
            dtype=object,
        )
    elif n == 10:
        p = np.array(
            [
                deltap[:9],
                deltap[9:17],
                deltap[17:24],
                deltap[24:30],
                deltap[30:35],
                deltap[35:39],
                deltap[39:42],
                deltap[42:44],
                deltap[44:],
            ],
            dtype=object,
        )
    else:
        log.critical("invalid hole number: %s", n)

    # calculates closure phases for general N-hole mask (with p-array set
    #     up properly above)
    cps = np.zeros((n - 1) * (n - 2) // 2)
    for j1 in range(n - 2):
        for j2 in range(n - 2 - j1):
            cps[int(j1 * ((n + (n - 3) - j1) / 2.0)) + j2] = (
                p[j1][0] + p[j1 + 1][j2] - p[j1][j2 + 1]
            )

    return cps


def closure_amplitudes(amps, n=7):
    """
    Calculate closure amplitudes.

    Parameters
    ----------
    amps : 1D float array
        Fringe amplitudes

    n : int, optional
        Number of holes (default=7)

    Returns
    -------
    CAs : 1D float array
        Closure amplitudes
    """
    arr = populate_symmamparray(amps, n=n)  # fringe amp array
    nn = 0

    cas = np.zeros(int(comb(n, 4)))

    for ii in range(n - 3):
        for jj in range(n - ii - 3):
            for kk in range(n - jj - ii - 3):
                for ll in range(n - jj - ii - kk - 3):
                    cas[nn + ll] = (
                        arr[ii, jj + ii + 1]
                        * arr[ll + ii + jj + kk + 3, kk + jj + ii + 2]
                        / (arr[ii, kk + ii + jj + 2] * arr[jj + ii + 1, ll + ii + jj + kk + 3])
                    )
                nn = nn + ll + 1

    return cas


def q4_phases(deltaps, n=7):
    """
    Calculate phases for each set of 4 holes.

    Parameters
    ----------
    deltaps : 1D float array
        Pistons between each pair of holes

    n : int, optional
        Number of holes (default=7)

    Returns
    -------
    quad_phases : 1D float array
        Quad phases
    """
    arr = populate_antisymmphasearray(deltaps, n=n)  # fringe phase array
    nn = 0
    quad_phases = np.zeros(int(comb(n, 4)))

    for ii in range(n - 3):
        for jj in range(n - ii - 3):
            for kk in range(n - jj - ii - 3):
                for ll in range(n - jj - ii - kk - 3):
                    quad_phases[nn + ll] = (
                        arr[ii, jj + ii + 1]
                        + arr[ll + ii + jj + kk + 3, kk + jj + ii + 2]
                        - arr[ii, kk + ii + jj + 2]
                        - arr[jj + ii + 1, ll + ii + jj + kk + 3]
                    )
                nn = nn + ll + 1

    return quad_phases
