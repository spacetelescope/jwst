import logging

import numpy as np

log = logging.getLogger(__name__)


class WeightedSigmaClip:
    """
    Iterative sigma-clipping with non-uniform inverse-variance weights.

    At each iteration the weighted mean and weighted standard deviation are
    computed, data points more than ``sigma_lower`` (below) or
    ``sigma_upper`` (above) standard deviations from the mean are masked, and
    the statistics are recomputed on the surviving points.  Iteration stops
    when the change in weighted mean falls below ``tolerance`` times the mean
    itself, or after ``iters`` iterations.

    Parameters
    ----------
    iters : int, optional
        Maximum number of clipping iterations.  Default is 10.
    tolerance : float, optional
        Convergence threshold.  Iteration stops when
        ``|Δmean| < tolerance * |mean|`` for all elements.  Default is 1e-3.
    sigma : float, optional
        Default clipping threshold in standard deviations, applied
        symmetrically when ``sigma_lower`` or ``sigma_upper`` are not
        specified.  Default is 3.0.
    sigma_lower : float or None, optional
        Lower clipping threshold in standard deviations (positive value;
        points more than this many sigma *below* the mean are clipped).
        If ``None``, defaults to ``sigma``.
    sigma_upper : float or None, optional
        Upper clipping threshold in standard deviations.  If ``None``,
        defaults to ``sigma``.
    """

    def __init__(self, iters=10, tolerance=1e-3, sigma=3.0, sigma_lower=None, sigma_upper=None):

        self.sigma_lower = -sigma_lower if sigma_lower is not None else -sigma
        self.sigma_upper = sigma_upper if sigma_upper is not None else sigma
        self.iters = iters
        self.tolerance = tolerance

    def __call__(self, data, weights, axis=None, return_mask=True):
        """
        Run weighted sigma-clipping on ``data``.

        Parameters
        ----------
        data : array_like
            Input data array.
        weights : array_like
            Non-negative inverse-variance weights, same shape as ``data``.
            NaN weights are treated as masked.  Negative values are clamped
            to zero.
        axis : int or None, optional
            Axis along which to clip.  If ``None`` (default), the statistic
            is computed over the entire flattened array.
        return_mask : bool, optional
            If ``True`` (default), return the final boolean mask as a third
            element.  Masked (clipped or originally invalid) pixels are
            ``True``.

        Returns
        -------
        ave : scalar or ndarray
            Weighted mean of the unclipped data.
        sig : scalar or ndarray
            Weighted standard deviation of the unclipped data.
        mask : ndarray of bool
            Boolean mask of clipped/invalid pixels (only returned when
            ``return_mask=True``).
        """
        # get new weights
        w = np.maximum(weights, 0.0)
        mask = np.isnan(data) | np.isnan(weights)

        # an internal variable
        x = np.ma.masked_array(data, mask=mask)

        # the initial setting for the output weighted mean
        ave = np.inf
        for _itr in range(self.iters):
            # a temporary mean for this iteration
            thisave = np.ma.average(x, weights=w, axis=axis, keepdims=True)

            # compute the difference, variance, and sigma
            dif = x - thisave
            var = np.ma.average(dif * dif, weights=w, axis=axis, keepdims=True)
            sig = np.sqrt(var)

            # scale differences by sigma
            z = dif / sig

            # update the mask
            x.mask |= (self.sigma_lower > z) | (z > self.sigma_upper)

            # check if we have converged
            dif = thisave - ave
            ave = thisave
            if np.all(np.abs(dif) < self.tolerance * np.abs(thisave)):
                break
        else:
            log.warning("Clipping did not converge.")

        ave = ave[0] if axis is None else np.squeeze(ave)
        sig = sig[0] if axis is None else np.squeeze(sig)

        if return_mask:
            return ave, sig, x.mask
        else:
            return ave, sig

    # def total(self, data, weights, **kwargs):
    #     """
    #     Return a sigma-clipped weighted sum.

    #     This is like

    #     .. math::
    #        total = sum(w*x)
    #     """
    #     kwargs["return_mask"] = True
    #     ave, _sig, _mask = self(data, weights, **kwargs)
    #     # tot = ave*np.nansum(mask*weights)
    #     return ave
