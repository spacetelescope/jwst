"""Fit a spectral shape to a WFSS dispersed source."""

import logging

import numpy as np

log = logging.getLogger(__name__)

__all__ = [
    "apply_basis_coeffs",
    "fit_slit_by_basis_images",
]


class SlitFitError(Exception):
    """Raise when spectral fitting fails."""

    pass


def _get_basis_images(simul_slit):
    """
    Collect the flat simulated image and all ``fluxmodel_N`` polynomial basis images.

    Parameters
    ----------
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit with attributes ``data`` (the constant/flat term) and optional
        ``fluxmodel_1``, ``fluxmodel_2``, ... attributes for higher-degree terms.

    Returns
    -------
    basis : list of ndarray
        ``[data, fluxmodel_1, fluxmodel_2, ...]`` as numpy arrays.
    """
    basis = [np.asarray(simul_slit.data)]
    k = 1
    while True:
        mc = getattr(simul_slit, f"fluxmodel_{k}", None)
        if mc is None:
            break
        basis.append(np.asarray(mc))
        k += 1
    return basis


def fit_slit_by_basis_images(observed_slit, simul_slit, l2_alpha=0.0, rejection_threshold=0.1):
    """
    Fit a linear combination of dispersed basis images to the observed slit.

    The constant (degree-0) term is ``simul_slit.data`` (the flat-spectrum simulation).
    Higher-degree terms are taken from the ``fluxmodel_1``, ``fluxmodel_2``, ...
    attributes of ``simul_slit``.  These are the grism-frame images produced by
    passing polynomial flux models through ``disperse()``.

    The fit solves::

        observed ≈ c_0 * data + c_1 * fluxmodel_1 + c_2 * fluxmodel_2 + ...

    via inverse-variance-weighted least squares on valid pixels, using the
    ``err`` array of the observed slit as pixel uncertainties.
    When ``l2_alpha > 0``, L2 regularisation is applied to the weighted normal equations.

    Parameters
    ----------
    observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Observed 2-D spectral cutout.
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit with ``data`` and optional ``fluxmodel_N`` attributes.
    l2_alpha : float, optional
        L2 regularisation strength.  Added to the diagonal of the weighted
        normal-equation matrix as ``alpha * I`` before solving, which penalizes
        large coefficients.  A value of ``0`` (the default) turns off regularization.
        Typical useful values are in the range ``1e-3`` - ``1e1``.
    rejection_threshold : float, optional
        If the fitted constant term coefficient ``c_0`` deviates from 1 by more than this amount,
        the fit is rejected and `None` is returned.  This fit rejection is necessary to
        avoid fits "blowing up" when a source is located in nonzero (pseudo-)background,
        either from a nearby bright source or because the background subtraction was imperfect.
        If None, no fits will be rejected.

    Returns
    -------
    coeffs : ndarray
        Best-fit coefficients ``[c_0, c_1, ...]``.

    Raises
    ------
    SlitFitError
        If there are fewer valid pixels than basis terms.
    """
    basis = _get_basis_images(simul_slit)
    obs_data = np.asarray(observed_slit.data)

    mask = np.isfinite(obs_data) & np.isfinite(basis[0]) & (basis[0] != 0)
    if getattr(observed_slit, "dq", None) is not None:
        mask &= (np.asarray(observed_slit.dq) & 1) == 0

    n_valid = int(mask.sum())
    n_terms = len(basis)
    if n_valid < n_terms:
        raise SlitFitError(
            f"Only {n_valid} valid pixel(s) available for a {n_terms}-term linear fit "
            f"(need at least {n_terms})."
        )

    # Build inverse-variance weights from the error array.
    err_arr = getattr(observed_slit, "err", None)
    if err_arr is not None:
        is_finite_err = np.isfinite(err_arr) & (err_arr > 0)
        if not np.any(is_finite_err & mask):
            raise SlitFitError(
                "No valid pixels have finite positive error values; cannot compute fit weights."
            )
        with np.errstate(divide="ignore", invalid="ignore"):
            inv_var = np.where(
                is_finite_err,
                1.0 / err_arr**2,
                0.0,
            )
        w_sqrt = np.sqrt(inv_var)[mask]
    else:
        w_sqrt = np.ones(n_valid)

    design_matrix = np.column_stack([b[mask] for b in basis])
    # Weighted normal equations: (A^T W A) c = A^T W b
    aw = design_matrix * w_sqrt[:, np.newaxis]
    bw = obs_data[mask] * w_sqrt
    if l2_alpha == 0.0:
        coeffs, *_ = np.linalg.lstsq(aw, bw, rcond=None)
    else:
        ata = aw.T @ aw
        atb = aw.T @ bw
        coeffs = np.linalg.solve(ata + l2_alpha * np.eye(n_terms), atb)

    # log some fit diagnostics for the source
    n_total = obs_data.size
    if rejection_threshold is not None and np.abs(coeffs[0] - 1) > rejection_threshold:
        log.debug(f"Fitted constant term c_0={coeffs[0]:.3g} is far from 1; rejecting fit.")
        return None
    log.debug(
        f"source_id={observed_slit.source_id} "
        f"order={observed_slit.meta.wcsinfo.spectral_order} "
        f"valid_pixels/total={n_valid}/{n_total} "
        # f"cond={cond:.3g} " # condition number of weighted design matrix
        f"coeffs={np.array2string(coeffs, precision=4, suppress_small=True)}"
    )
    return coeffs


def apply_basis_coeffs(simul_slit, coeffs):
    """
    Reconstruct a fitted slit as a linear combination of dispersed basis images.

    Parameters
    ----------
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit with ``data`` and optional ``fluxmodel_N`` attributes.
    coeffs : ndarray
        Coefficients from `fit_slit_by_basis_images`.

    Returns
    -------
    ndarray
        Fitted slit image.
    """
    basis = _get_basis_images(simul_slit)
    return sum(c * b for c, b in zip(coeffs, basis, strict=True))
