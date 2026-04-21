"""Fit a spectral shape to a WFSS slit."""

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
        Simulated slit carrying ``data`` (the constant/flat term) and optional
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


def fit_slit_by_basis_images(observed_slit, simul_slit):
    """
    Fit a linear combination of dispersed basis images to the observed slit.

    The constant (degree-0) term is ``simul_slit.data`` (the flat-spectrum simulation);
    higher-degree terms are taken from the ``fluxmodel_1``, ``fluxmodel_2``, ...
    attributes on ``simul_slit``.  These are the grism-frame images produced by
    passing polynomial flux models (``λ``, ``λ²``, ...) through ``disperse()``.

    The fit solves::

        observed ≈ c_0 * data + c_1 * fluxmodel_1 + c_2 * fluxmodel_2 + ...

    via linear least-squares on valid pixels.

    Parameters
    ----------
    observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Calibrated observed 2-D spectral cutout.
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit with ``data`` and optional ``fluxmodel_N`` attributes.

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

    design_matrix = np.column_stack([b[mask] for b in basis])
    coeffs, *_ = np.linalg.lstsq(design_matrix, obs_data[mask], rcond=None)
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
