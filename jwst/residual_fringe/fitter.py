import numpy as np

import scipy.interpolate


def _lsq_spline(x, y, weights, knots, degree):
    return scipy.interpolate.LSQUnivariateSpline(x, y, knots, w=weights, k=degree)


def spline_fitter(x, y, weights, knots, degree, reject_outliers=False, domain=10, tolerance=0.0001):
    if not reject_outliers:
        return _lsq_spline(x, y, weights, knots, degree)

    # fit with chi sq outlier rejection
    # helpers
    def chi_sq(spline, weights):
        return np.nansum((y - spline(x)) ** 2 * weights)


    # initial fit
    spline = _lsq_spline(x, y, weights, knots, degree)
    chi = chi_sq(spline, weights)

    # astropy code used the model params which pad the knots based on degree
    nparams = len(knots) + (degree + 1) * 2
    deg_of_freedom = np.sum(weights) - nparams

    for _ in range(1000 * nparams):
        scale = np.sqrt(chi / deg_of_freedom)

        # Calculate new weights
        resid = (y - spline(x)) / (scale * domain)
        new_w = (np.where(resid**2 <= 1, 1 - resid**2, 0.))**2 * weights

        # Fit new model and find chi
        spline = _lsq_spline(x, y, new_w, knots, degree)
        new_chi = chi_sq(spline, new_w)

        # Check if fit has converged
        tol = tolerance if new_chi < 1 else tolerance * new_chi
        if np.abs(chi - new_chi) < tol:
            break
        chi = new_chi
    else:
        raise RuntimeError("Bad fit, method should have converged")

    return spline
