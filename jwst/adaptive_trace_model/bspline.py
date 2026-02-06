import logging

import numpy as np
from astropy.stats import sigma_clipped_stats as scs
from scipy.interpolate import make_lsq_spline

__all__ = ["bspline_fit"]

log = logging.getLogger(__name__)


def bspline_fit(
    xvec,
    yvec,
    nbkpts=50,
    wrapsig_low=3.0,
    wrapsig_high=3.0,
    wrapiter=3,
    space_ratio=1.2,
    verbose=False,
):
    """
    Fit a univariate basis spline to a vector with iterative rejection.

    Parameters
    ----------
    xvec : ndarray
        1D array of x values for the vector.
    yvec : ndarray
        1D array of y values to fit.
    nbkpts : int, optional
        Number of spline breakpoints (knots).
    wrapsig_low : float, optional
        Low sigma threshold for iterative fit.
    wrapsig_high : float, optional
        High sigma threshold for iterative fit.
    wrapiter : int, optional
        Number of iterations for the fit.
    space_ratio : float, optional
        Maximum spacing ratio to allow fitting to continue. If
        the tenth-largest spacing in the input ``xvec`` is larger
        than the knot spacing by this ratio, then return None instead
        of attempting to fit the data.
    verbose : bool, optional
        If True, an info log message is generated for each fitting iteration.

    Returns
    -------
    spline_model : `~scipy.interpolate.BSpline` or None
        The spline model fit to the data.  If the fit failed, None is returned.
    """
    spline_model = None
    indx = np.isfinite(xvec) & np.isfinite(yvec)
    xvec_use = xvec[indx]
    yvec_use = yvec[indx]

    # Reject points with abnormally large spacing as they can
    # cause problems with breakpoints
    spacing = np.abs(np.diff(xvec_use, prepend=0))
    space_mean, _, space_rms = scs(spacing)
    if np.isfinite(space_rms):
        good = spacing <= (space_mean + 5 * space_rms)
        xvec_use = xvec_use[good]
        yvec_use = yvec_use[good]

    # If the tenth-largest spacing was bigger than the knot spacing
    # (with some margin) then don't bspline
    knotspacing = (np.max(xvec_use) - np.min(xvec_use)) / nbkpts
    tenthspace = np.partition(spacing, -10)[-10]
    if tenthspace > (space_ratio * knotspacing):
        return spline_model

    # Number of points before iterative loop
    norig = len(xvec_use)

    last_spline = None
    for ii in range(0, wrapiter):
        knotspacing = (np.max(xvec_use) - np.min(xvec_use)) / nbkpts
        knotmin = np.min(xvec_use) + knotspacing / 2.0
        knotmax = np.max(xvec_use) - knotspacing / 2.0
        knots = np.arange(knotmin, knotmax, knotspacing)

        # Add exterior knots for cubic spline
        degree = 3
        xb = xvec_use[0]
        xe = xvec_use[-1]
        knots = np.concatenate(([xb] * (degree + 1), knots, [xe] * (degree + 1)))

        try:
            spline_model = make_lsq_spline(xvec_use, yvec_use, knots, k=degree)
        except ValueError:
            # Bad fit: break the loop and return the last fit spline (or None)
            spline_model = last_spline
            break
        ytemp = spline_model(xvec_use)
        if np.any(np.isnan(ytemp)):
            # Bad fit: break the loop and return the last fit spline (or None)
            spline_model = last_spline
            break

        # Good fit: keep it for the next loop.
        last_spline = spline_model

        # No need to reject data on the last iteration
        if ii == wrapiter - 1:
            break

        # If continuing on, reject significant outliers
        diff = yvec_use - ytemp
        rms = np.std(diff)
        rej = (diff < -wrapsig_low * rms) | (diff > wrapsig_high * rms)
        keep = ~rej
        if verbose:
            log.info(f"Iter {ii} Rejected {np.sum(rej)} Kept {np.sum(keep)}")

        # If no points rejected or too few kept break out of the loop
        if np.sum(keep) < 0.8 * norig:
            # Too many rejected - keep the last model instead
            if verbose:
                log.info("Stopping iteration: too many points rejected")
            spline_model = last_spline
            break
        elif np.sum(rej) == 0:
            break

        # Otherwise update the x and y vectors and try again
        xvec_use = xvec_use[keep]
        yvec_use = yvec_use[keep]

    return spline_model
