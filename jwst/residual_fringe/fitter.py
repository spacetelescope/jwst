import numpy as np

from astropy.modeling.fitting import _model_to_fit_params, LevMarLSQFitter


class DegenerateEvidenceError(Exception):
    """
    Raised when evidence calculation becomes degenerate.
    """


class ChiSqOutlierRejectionFitter:
    """Chi squared statistic for outlier rejection"""

    def __init__(self, fitter, domain=None, tolerance=0.0001):
        self.fitter = fitter
        self.tolerance = tolerance

        if domain is None:
            self.domain = 10
        else:
            self.domain = domain

    @staticmethod
    def kernel(x, weights=None):
        """
        Weighting function dependent only on provided value (usualy a residual)
        """

        kernal = (np.where(x**2 <= 1, 1 - x**2, 0.))**2
        if weights is not None:
            kernal *= weights

        return kernal

    @staticmethod
    def _chi(model, x, y, weights=None):

        resid = (y - model(x))**2
        if weights is not None:
            resid *= weights

        return np.sum(resid)

    @staticmethod
    def _params(model):
        return _model_to_fit_params(model)[0]

    @staticmethod
    def _sum_weights(x, weights=None):
        if weights is None:
            return len(x)
        else:
            return np.sum(weights)

    def _deg_of_freedom(self, model, x, weights=None):
        nparams = len(self._params(model))
        sum_weights = self._sum_weights(x, weights)

        return nparams, sum_weights - nparams

    @staticmethod
    def _scale(chi, deg_of_freedom):
        return np.sqrt(chi / deg_of_freedom)

    def _log_likelihood(self, chi, x, weights=None, var=1.0):
        sum_weights = self._sum_weights(x, weights)

        return -0.5 * (sum_weights * np.log(2 * np.pi * var) + chi / var)

    def _make_variance(self, chi, deg_of_freedom, scale=None):
        if scale is None:
            scale = self._scale(chi, deg_of_freedom)

        return scale ** 2

    @staticmethod
    def _jacobian(model, x, weights=None):
        jacobian = np.array(model.jacobian(x, *model.parameters))
        if weights is not None:
            jacobian = np.array([np.ravel(this)
                                 for this in np.array(weights) * jacobian])

        return jacobian

    def _hessian(self, model, x, weights=None):

        jacobian = self._jacobian(model, x, weights)

        if weights is None:
            return 2 * np.inner(jacobian, jacobian)
        else:
            return 2 * np.inner(jacobian, jacobian * weights)

    def _log(self, value):
        err_state = np.seterr(all='raise')
        try:
            log = np.log(value)
        except FloatingPointError:
            np.seterr(**err_state)

            raise DegenerateEvidenceError("The evidence calculation has become degenerate")

        np.seterr(**err_state)

        return log

    def _get_log_z(self, model, x, y, weights,
                   limits, noise_limits, fixed_scale):
        chi = self._chi(model, x, y, weights)
        nparams, deg_of_freedom = self._deg_of_freedom(model, x, weights)

        if deg_of_freedom <= 0:
            raise RuntimeError("More parameters than (weighted) data points")

        if noise_limits is None:
            scale = 1.0 if fixed_scale is None else fixed_scale
            var = scale**2
            log_occam = 0.0
            spr = 0.0
        else:
            scale_range = np.log(noise_limits[1]) - np.log(noise_limits[0])
            var = self._make_variance(chi, deg_of_freedom)
            log_occam = 0.5 * np.log(np.pi * deg_of_freedom) - scale_range
            spr = np.log(scale_range)

        log_likelihood = self._log_likelihood(chi, x, weights, var)

        if nparams == 0:
            return log_likelihood + log_occam

        prior_range = np.atleast_1d(np.log(np.array(limits[1]) - np.array(limits[0])))

        prior_length = len(prior_range)
        if prior_length == 1:
            spr += prior_range[0]
        else:
            spr += np.sum(prior_range)

        if prior_length < nparams:
            spr += (nparams - prior_length) * prior_range[-1]

        det_hessian = np.linalg.det(self._hessian(model, x, weights))
        log_det = self._log(det_hessian)

        log_occam += -spr + 0.5 * (nparams * self._log(2 * np.pi * var) - log_det)

        log_z = log_likelihood + log_occam

        return log_z

    def _get_evidence(self, model, x, y, weights,
                      limits, noise_limits, fixed_scale):
        return self._get_log_z(model, x, y, weights,
                               limits, noise_limits, fixed_scale) / np.log(10)

    def __call__(self, model, x, y, weights=None, get_evidence=False, **kwargs):
        # Assume equal weights if none are provided

        limits = kwargs.pop('limits')
        noise_limits = kwargs.pop('noise_limits', None)
        fixed_scale = kwargs.pop('fixed_scale', None)

        new_model = model.copy()

        # perform the initial fit
        new_model = self.fitter(new_model, x, y, weights=weights, **kwargs)
        chi = self._chi(new_model, x, y, weights)

        # calculate degrees of freedom
        nparams, deg_of_freedom = self._deg_of_freedom(new_model, x, weights)

        # Iteratively adjust the weights until fit converges
        for _ in range(1000 * nparams):
            scale = self._scale(chi, deg_of_freedom)

            # Calculate new weights
            resid = (y - new_model(x)) / (scale * self.domain)
            new_weights = self.kernel(resid, weights)

            # Fit new model and find chi
            new_model = self.fitter(new_model, x, y, weights=new_weights, **kwargs)
            new_chi = self._chi(new_model, x, y, new_weights)

            # Check if fit has converged
            tol = self.tolerance if new_chi < 1 else self.tolerance * new_chi
            if np.abs(chi - new_chi) < tol:
                break
            chi = new_chi
        else:
            raise RuntimeError("Bad fit, method should have converged")

        if get_evidence:
            evidence = self._get_evidence(new_model, x, y, new_weights,
                                          limits, noise_limits, fixed_scale)

            return new_model, evidence, chi
        else:
            return new_model
