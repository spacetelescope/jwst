import numpy as np

from astropy.modeling.fitting import model_to_fit_params


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

        return np.nansum(resid)

    @staticmethod
    def _params(model):
        return model_to_fit_params(model)[0]

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

    def __call__(self, model, x, y, weights=None, **kwargs):
        # Assume equal weights if none are provided

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
            new_w = self.kernel(resid, weights)

            # Fit new model and find chi
            new_model = self.fitter(new_model, x, y, weights=new_w, **kwargs)
            new_chi = self._chi(new_model, x, y, new_w)

            # Check if fit has converged
            tol = self.tolerance if new_chi < 1 else self.tolerance * new_chi
            if np.abs(chi - new_chi) < tol:
                break
            chi = new_chi
        else:
            raise RuntimeError("Bad fit, method should have converged")

        return new_model
