"""Fit a spectral shape to a WFSS slit."""

import logging

import matplotlib.pyplot as plt
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


def fit_slit_by_basis_images(observed_slit, simul_slit, l2_alpha=0.0):
    """
    Fit a linear combination of dispersed basis images to the observed slit.

    The constant (degree-0) term is ``simul_slit.data`` (the flat-spectrum simulation);
    higher-degree terms are taken from the ``fluxmodel_1``, ``fluxmodel_2``, ...
    attributes on ``simul_slit``.  These are the grism-frame images produced by
    passing polynomial flux models (``λ``, ``λ²``, ...) through ``disperse()``.

    The fit solves::

        observed ≈ c_0 * data + c_1 * fluxmodel_1 + c_2 * fluxmodel_2 + ...

    via ridge regression (L2-regularised least squares) on valid pixels.
    When ``l2_alpha=0`` (the default) this reduces to ordinary least squares.

    Parameters
    ----------
    observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Calibrated observed 2-D spectral cutout.
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit with ``data`` and optional ``fluxmodel_N`` attributes.
    l2_alpha : float, optional
        L2 regularisation strength.  Added to the diagonal of the normal-equation
        matrix as ``alpha * I`` before solving, which penalises large coefficients
        and suppresses the oscillating sign-alternating solutions that arise when
        the monomial basis images are nearly collinear.  A value of ``0`` (the
        default) gives ordinary least squares.  Typical useful values are in the
        range ``1e-3`` – ``1e1``; the optimal value is data-dependent.

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
    if l2_alpha == 0.0:
        coeffs, *_ = np.linalg.lstsq(design_matrix, obs_data[mask], rcond=None)
    else:
        ata = design_matrix.T @ design_matrix
        atb = design_matrix.T @ obs_data[mask]
        coeffs = np.linalg.solve(ata + l2_alpha * np.eye(n_terms), atb)
    # _basis_diagnostic_plot(observed_slit, basis)
    # _diagnostic_plot(observed_slit, simul_slit, basis, coeffs)
    return coeffs


def _basis_diagnostic_plot(observed_slit, basis):
    """
    Show the ratio of each higher-order basis image to the 0th-order basis.

    Layout: 1 row × (n_basis - 1) columns, all panels on a shared color scale.
    - Each panel: basis[k] / basis[0], for k = 1 .. n_basis-1
    """
    n_basis = len(basis)
    if n_basis < 2:
        return  # nothing to compare

    sid = getattr(observed_slit, "source_id", "?")
    order = getattr(
        getattr(getattr(observed_slit, "meta", None), "wcsinfo", None),
        "spectral_order",
        "?",
    )

    ref = np.asarray(basis[0], dtype=float)
    ref_safe = np.where(ref == 0, np.nan, ref)  # avoid div-by-zero

    # Compute all ratios first so we can find a shared color scale
    ratios = [np.asarray(basis[k], dtype=float) / ref_safe for k in range(1, n_basis)]
    all_deviations = np.concatenate([np.abs(r - 1)[np.isfinite(r)] for r in ratios])
    rlim = np.percentile(all_deviations, 99) if len(all_deviations) else 1.0

    n_cols = n_basis - 1
    fig, axes = plt.subplots(1, n_cols, figsize=(3 * n_cols, 12), squeeze=False)
    fig.suptitle(f"basis ratios — src={sid} order={order}", fontsize=11)

    last_im = None
    for col_idx, ratio in enumerate(ratios):
        k = col_idx + 1
        last_im = axes[0, col_idx].imshow(
            ratio, origin="lower", vmin=1 - rlim, vmax=1 + rlim, aspect="auto", cmap="RdBu_r"
        )
        axes[0, col_idx].set_title(f"basis {k} / basis 0")
        if col_idx > 0:
            axes[0, col_idx].tick_params(labelleft=False)

    fig.colorbar(last_im, ax=axes[0, :], location="bottom", orientation="horizontal", shrink=0.6)

    plt.tight_layout()
    plt.show()


def _diagnostic_plot(observed_slit, _simul_slit, basis, coeffs):
    """
    Show a diagnostic figure for one polynomial basis fit.

    Panels (top to bottom):
    - observed slit
    - each basis term (degree 0, 1, 2, ...) scaled by its fitted coefficient
    - combined fitted simulation
    - residual (observed - fitted)
    """
    sid = getattr(observed_slit, "source_id", "?")
    order = getattr(
        getattr(getattr(observed_slit, "meta", None), "wcsinfo", None),
        "spectral_order",
        "?",
    )
    # print(
    #     f"[wavefit] source_id={sid} order={order} "
    #     f"coeffs=[{', '.join(f'{c:.4g}' for c in coeffs)}]"
    # )

    obs_data = np.asarray(observed_slit.data)
    fitted = sum(c * b for c, b in zip(coeffs, basis, strict=True))
    residual = obs_data - fitted

    n_basis = len(basis)
    n_panels = 1 + n_basis + 1 + 1  # obs + each basis term + fitted + residual
    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 8), squeeze=False)
    fig.suptitle(f"wavefit diagnostic — src={sid} order={order}", fontsize=11)

    vmin = np.nanpercentile(obs_data, 1)
    vmax = np.nanpercentile(obs_data, 99)
    dlim = np.nanpercentile(np.abs(residual), 99)

    ax_idx = 0

    # Panel 0: observed slit
    axes[0, ax_idx].imshow(
        obs_data, origin="lower", vmin=vmin, vmax=vmax, aspect="auto", cmap="gray_r"
    )
    axes[0, ax_idx].set_title("observed")
    fig.colorbar(
        axes[0, ax_idx].images[0], ax=axes[0, ax_idx], location="bottom", orientation="horizontal"
    )
    ax_idx += 1

    # Panels 1..n_basis: each scaled basis term
    for k, (c, b) in enumerate(zip(coeffs, basis, strict=True)):
        scaled = c * b
        bvmin = np.nanpercentile(scaled, 1)
        bvmax = np.nanpercentile(scaled, 99)
        axes[0, ax_idx].imshow(
            scaled, origin="lower", vmin=bvmin, vmax=bvmax, aspect="auto", cmap="gray_r"
        )
        axes[0, ax_idx].set_title(f"basis term {k} × coeff {c:.4g}")
        axes[0, ax_idx].tick_params(labelleft=False)
        fig.colorbar(
            axes[0, ax_idx].images[0],
            ax=axes[0, ax_idx],
            location="bottom",
            orientation="horizontal",
        )
        ax_idx += 1

    # Panel: combined fitted simulation
    axes[0, ax_idx].imshow(
        fitted, origin="lower", vmin=vmin, vmax=vmax, aspect="auto", cmap="gray_r"
    )
    axes[0, ax_idx].set_title("fitted simulation (combined)")
    axes[0, ax_idx].tick_params(labelleft=False)
    fig.colorbar(
        axes[0, ax_idx].images[0], ax=axes[0, ax_idx], location="bottom", orientation="horizontal"
    )
    ax_idx += 1

    # Panel: residual (observed - fitted)
    im_res = axes[0, ax_idx].imshow(
        residual, origin="lower", vmin=-dlim, vmax=dlim, aspect="auto", cmap="RdBu_r"
    )
    axes[0, ax_idx].set_title(f"residual (obs - fitted)  std={np.nanstd(residual):.4g}")
    axes[0, ax_idx].tick_params(labelleft=False)
    fig.colorbar(im_res, ax=axes[0, ax_idx], location="bottom", orientation="horizontal")

    plt.tight_layout()
    plt.show()


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
