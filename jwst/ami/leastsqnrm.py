import logging
import warnings

import numpy as np
import numpy.linalg as linalg
from scipy.special import comb


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


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
        Bad pixel mask, same shape as image.

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

    if dqm is not None:
        flatdqm = dqm.reshape(np.shape(img)[0] * np.shape(img)[1])
        nanlist = np.where(flatdqm)  # where DO_NOT_USE up.
    else:
        nanlist = (np.array((), dtype=np.int32),)  # shouldn't occur w/MAST JWST data

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


def matrix_operations(img, model, linfit=False, dqm=None):
    """
    Use least squares matrix operations to solve A x = b.

    A is the model, b is the data (img), and x is the coefficient vector we are solving for.
    In 2-D, data x = inv(At.A).(At.b).

    TODO: replace linearfit with scipy fitting.

    Parameters
    ----------
    img : 2D float array
        Input data
    model : 2D float array
        Analytic model
    linfit : bool
        Whether to perform linear fit
    dqm : 2D bool array
        Bad pixel mask slice, same shape as image.

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
    log.info("fringefitting.leastsqnrm.matrix_operations(): ")
    log.info(f"\timg {img.shape:}")
    log.info(
        f"\tL x W = {img.shape[0]:d} x {img.shape[1]:d} = {img.shape[0] * img.shape[1]:d}",
    )
    log.info(f"\tflatimg {flatimg.shape:}")

    log.info("")
    log.info("\ttype(dqm) %s", type(dqm))
    if dqm is not None:
        flatdqm = dqm.reshape(np.shape(img)[0] * np.shape(img)[1])
        nanlist = np.where(flatdqm)  # where DO_NOT_USE up.
    else:
        nanlist = (np.array((), dtype=np.int32),)  # shouldn't occur w/MAST JWST data

    log.info(f"\ttype(nanlist) {type(nanlist):}, len={len(nanlist):}")
    log.info(f"\tnumber of nanlist pixels: {len(nanlist[0]):d} items")
    log.info(f"\t{len(nanlist[0]):d} DO_NOT_USE pixels found in data slice")

    flatimg = np.delete(flatimg, nanlist)

    log.info(f"\tflatimg {flatimg.shape:} after deleting {len(nanlist[0]):d}")

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

    log.info("data flux %s", flatimg.sum())
    log.info("flat model dimensions %s", np.shape(flatmodel))
    log.info("model transpose dimensions %s", np.shape(flatmodeltransp))
    log.info("flat image dimensions %s", np.shape(flatimg))
    log.info("transpose * image data dimensions %s", np.shape(data_vector))
    log.info("flat img * transpose dimensions %s", np.shape(inverse))

    if linfit:
        # photon noise
        noise = np.sqrt(np.abs(flatimg))

        # this sets the weights of pixels fulfilling condition to zero
        weights = np.where(np.abs(flatimg) <= 1.0, 0.0, 1.0 / (noise**2))

        # initialize object with uniform weights
        result = LinearFit(flatimg, np.diag(weights), flatmodeltransp)

        # do the fit
        result.fit()

        # delete inverse_covariance_matrix to reduce size of pickled file
        result.inverse_covariance_matrix = []

        linfit_result = result
        log.info("Returned linearfit result")

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


class LinearFit:
    """
    Perform a general least-squares fit of a linear model using numpy matrix inversion.

    Uncertainties in the dependent variables (but not in the independent variables)
    can be taken into account.
    All inputs have to be numpy matrices.

    Math is based on Press'
    Numerical Recipes p661 : Section 15.2 Fitting Data to a Straight Line
    Numerical Recipes p671 : Section 15.4 General Linear Least Squares

    Code is based on an early yorick implementation by Damien Segransan (University of Geneva)
    Python implementation and tools by Johannes Sahlmann 2009-2017
    (University of Geneva, European Space Agency, STScI/AURA)

    Attributes
    ----------
    dependent_variable : np.ndarray (1xN)
        Dependent_variables of the linear equation system (N equations, M unknown coefficients)
    inverse_covariance_matrix : np.ndarray (NxN)
        Inverse covariance matrix corresponding to the dependent_variable.
        i.e. data weights proportional to 1/sigma**2 where sigma=uncertainty
    independent_variable : np.ndarray (MxN)
        The independent_variables that are multiplied by the unknown coefficients

    Calculated Attributes
    ----------
    p : np.ndarray
        Coefficients of the solution
    p_formal_uncertainty : np.ndarray
        Formal uncertainty of the coefficients
    p_formal_covariance_matrix : np.ndarray
        Formal covariance matrix of the coefficients (not rescaled)
    p_normalised_uncertainty : np.ndarray
        Normalised uncertainty (chi2 = 1) of the coefficients
    p_normalised_covariance_matrix : np.ndarray
        Normalised covariance matrix of the coefficients (rescaled to yield chi2=1)
    p_correlation_matrix : np.ndarray
        Coefficient correlation matrix
    fit : np.ndarray
        Values of the best-fit model
    residuals : np.ndarray
        Observed - Calculated (O-C) residuals
    chi2 : float
        Chi-square value of the best fit
    """

    def __init__(self, dependent_variable, inverse_covariance_matrix, independent_variable):
        self.y_i = dependent_variable
        self.X_ij = independent_variable
        self.inverse_covariance_matrix = inverse_covariance_matrix

    def fit(self):
        """Perform the linear fit."""
        # number of measurements/constraints, i.e. number of equations
        self.n_measurement = self.X_ij.shape[1]
        # number of coefficients, i.e. free parameters
        self.n_param = self.X_ij.shape[0]
        # number of degrees of freedom
        self.n_freedom = self.n_measurement - self.n_param

        a_ij = (self.X_ij) @ self.inverse_covariance_matrix
        alpha_kj = a_ij @ (self.X_ij.T)
        beta_k = a_ij @ self.y_i.T

        a_j = np.linalg.solve(alpha_kj, beta_k)

        yfit_i = (self.X_ij.T * a_j).T
        chi2 = ((yfit_i - self.y_i) @ self.inverse_covariance_matrix) @ (yfit_i - self.y_i).T
        c_jk = np.linalg.inv(alpha_kj)
        var_aj = np.diag(c_jk)

        # the variance on the fitted coefficients are the diagonal terms
        # of the coefficient covariance matrix
        # if the error bars are not well estimated, it is useful to
        # renormalize the variance by the measured chi2
        # divided by the expected chi2.
        # That means the normalised variances take into account the residual dispersion in the data
        normalized_variance_aj = var_aj.T * chi2 / self.n_freedom
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=RuntimeWarning, message="invalid value encountered in sqrt"
            )
            stdev_aj = np.sqrt(normalized_variance_aj)

        # coefficients of the solution
        self.p = np.array(a_j).flatten()

        # normalised uncertainty (chi2 = 1) of the coefficients
        self.p_normalised_uncertainty = np.array(stdev_aj).flatten()

        # values of the best-fit model
        self.fit_values = np.array(yfit_i).flatten()

        # Observed - Calculated (O-C) residuals
        self.residuals = np.array(self.y_i - yfit_i).flatten()

        # chi-square value of the best fit
        self.chi2 = np.array(chi2)[0][0]

        # formal uncertainty  of the coefficients
        self.p_formal_uncertainty = np.array(np.sqrt(var_aj)).flatten()

        # formal covariance matrix of the coefficients (not rescaled)
        self.p_formal_covariance_matrix = c_jk

        # normalised covariance matrix of the coefficients (rescaled to yield chi2=1)
        self.p_normalised_covariance_matrix = c_jk * self.chi2 / self.n_freedom

        #     compute correlation Matrix
        tmp_v = 1.0 / self.p_formal_uncertainty
        tmp = np.vstack((tmp_v, np.tile(np.zeros(len(tmp_v)), (len(tmp_v) - 1, 1))))
        m = tmp.T
        err_matrix = m.dot(m.T)
        correlation_matrix = np.multiply(err_matrix, self.p_formal_covariance_matrix.T)
        self.p_correlation_matrix = correlation_matrix
