"""Fit a spectral shape to a WFSS slit."""

import logging

import numpy as np

log = logging.getLogger(__name__)

__all__ = ["SlitIterativePolynomialFitter", "SlitPolynomialFitter", "apply_flam_to_slit"]


class SlitFitError(Exception):
    """Raise when spectral fitting fails."""

    pass


def _build_fit_arrays(observed_slit, simul_slit):
    """
    Extract masked pixel arrays needed for polynomial fitting.

    Parameters
    ----------
    observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        The calibrated observed 2-D spectral cutout.
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        The flat-spectrum simulated cutout. Must carry a ``wavelength`` array
        with the same shape as ``data``.

    Returns
    -------
    data_masked : ndarray, shape (N,)
        Observed pixel values at valid locations.
    sim_masked : ndarray, shape (N,)
        Simulated pixel values at valid locations.
    lam : ndarray, shape (N,)
        Wavelength values at valid locations.
    lam_ref : float
        Median wavelength, used as the polynomial centering point.

    Raises
    ------
    SlitFitError
        If ``simul_slit.wavelength`` is ``None`` or has the wrong shape.
    """
    data = np.asarray(observed_slit.data)
    sim = np.asarray(simul_slit.data)
    wavelength = np.asarray(simul_slit.wavelength)

    if wavelength is None or wavelength.shape != sim.shape:
        raise SlitFitError(
            "simul_slit.wavelength must be a 2-D array with the same shape as simul_slit.data"
        )

    dq = np.asarray(observed_slit.dq, dtype=np.uint32)
    mask = (
        (sim != 0)  # simulation has signal here
        & (wavelength > 0)  # wavelength was assigned
        & np.isfinite(data)  # observed pixel is valid
        & np.isfinite(sim)  # simulated pixel is valid
        & ((dq & 1) == 0)  # dq bit 0 is DO_NOT_USE
    )

    if not np.any(mask):
        raise SlitFitError(
            "No valid pixels found after masking. Check that simul_slit has nonzero signal, "
            "valid wavelengths, and that observed_slit has finite data with no DO_NOT_USE pixels."
        )

    lam_ref = float(np.median(wavelength[mask]))
    return data[mask], sim[mask], wavelength[mask], lam_ref


def _make_poly(coeffs, lam_ref):
    """
    Build a callable polynomial from lstsq coefficients.

    Parameters
    ----------
    coeffs : ndarray
        Polynomial coefficients ``[c_0, c_1, ..., c_d]`` such that
        ``p(λ) = Σ_k c_k * (λ - lam_ref)^k``.
    lam_ref : float
        Wavelength centering constant.

    Returns
    -------
    callable
        A function ``f(wavelength)`` that evaluates the polynomial.
    """

    def f(wavelength):
        wavelength = np.asarray(wavelength)
        d = wavelength - lam_ref
        poly_vals = np.zeros(wavelength.shape, dtype=float)
        for k, c in enumerate(coeffs):
            poly_vals += c * d**k
        return poly_vals

    return f


class SlitPolynomialFitter:
    """
    Fit a single polynomial spectral shape to the slit.

    Parameters
    ----------
    degree : int, optional
        Degree of the fitting polynomial.  Default is 2.
    """

    def __init__(self, degree=2):
        self.degree = degree

    def __call__(self, observed_slit, simul_slit):
        """
        Scale the simulation by a best-fit polynomial spectral shape.

        The flat-spectrum ``simul_slit.data`` is scaled pixel-by-pixel by a
        polynomial in wavelength::

            scaled[y, x] = simul[y, x] * p(λ[y, x])
            p(λ) = Σ_k  c_k * (λ - λ_ref)^k

        The coefficients are chosen via linear least squares to minimise the
        squared residuals between the scaled simulation and the observed data::

            min  ||observed.data  -  simul.data * p(simul.wavelength)||²

        Because ``p`` is linear in the coefficients, the design matrix columns
        are ``simul.data * (λ - λ_ref)^k``, so no iterative solver is needed.

        Parameters
        ----------
        observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
            The calibrated observed 2-D spectral cutout.
        simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
            The flat-spectrum simulated cutout for the same source and order.
            Must have a ``wavelength`` array with the same shape as ``data``.

        Returns
        -------
        f_lam : callable
            A function ``f_lam(wavelength)`` that evaluates the best-fit
            polynomial p(λ).

        Raises
        ------
        SlitFitError
            If ``simul_slit.wavelength`` is ``None`` or has a different shape
            from ``simul_slit.data``.
        SlitFitError
            If the number of valid pixels (after masking) is smaller than
            ``self.degree + 1``.

        Notes
        -----
        A pixel is included in the fit only when all of the following hold:

        * ``simul_slit.data != 0``  (simulation has signal here)
        * ``simul_slit.wavelength > 0``  (wavelength was assigned)
        * ``np.isfinite(observed_slit.data)``  (observed pixel is valid)
        * ``np.isfinite(simul_slit.data)``  (simulated pixel is valid)
        * DQ bit 0 (DO_NOT_USE) is not set in ``observed_slit.dq``
        """
        degree = self.degree
        data_masked, sim_masked, lam, lam_ref = _build_fit_arrays(observed_slit, simul_slit)

        n_valid = len(data_masked)
        if n_valid < degree + 1:
            raise SlitFitError(
                f"Only {n_valid} valid pixel(s) available for a degree-{degree} polynomial fit "
                f"(need at least {degree + 1}). Reduce the fitting degree or check the input data."
            )

        dlam = lam - lam_ref
        design_matrix = np.column_stack([sim_masked * dlam**k for k in range(degree + 1)])
        coeffs, *_ = np.linalg.lstsq(design_matrix, data_masked, rcond=None)
        return _make_poly(coeffs, lam_ref)


class SlitIterativePolynomialFitter:
    """
    Fit a polynomial spectral shape to a WFSS slit by incrementally adding terms.

    Starting from a constant (degree-0) scaling, one coefficient at a time is added
    up to ``max_degree``.  At each step ``k`` the previously fitted coefficients are
    held fixed and only the new ``c_k * (λ - λ_ref)^k`` term is determined by a
    1-D linear least-squares solve::

        c_k = (sim * (λ - λ_ref)^k) · residual_k
              ─────────────────────────────────────
              ‖ sim * (λ - λ_ref)^k ‖²

    where ``residual_k = data - sim * p_{k-1}(λ)``.  The result is a single
    polynomial of degree at most ``max_degree``.

    Because each step is a 1-D solve, the previous solution acts as a warm start:
    each new term explains only the multiplicative residual left by all lower-degree
    terms, without re-fitting those coefficients.

    Parameters
    ----------
    max_degree : int
        Maximum degree of the fitted polynomial.
    improvement_threshold : float or None, optional
        Minimum required relative improvement in the residual sum-of-squares
        (RSS) for each successive term (degree ≥ 1) to be accepted.  Computed as
        ``(RSS_before - RSS_after) / RSS_before``.  The constant term (degree 0)
        is always fitted.  If a higher-degree term does not meet the threshold,
        iteration stops and the polynomial built so far is returned.  ``None``
        (default) disables early stopping.
    """

    def __init__(self, max_degree, improvement_threshold=None):
        self.max_degree = max_degree
        self.improvement_threshold = improvement_threshold

    def __call__(self, observed_slit, simul_slit):
        """
        Fit a polynomial spectral shape by incrementally adding one term at a time.

        Parameters
        ----------
        observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
            The calibrated observed 2-D spectral cutout.
        simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
            The flat-spectrum simulated cutout for the same source and order.
            Must have a ``wavelength`` array with the same shape as ``data``.

        Returns
        -------
        f_lam : callable
            A function ``f_lam(wavelength)`` that evaluates the best-fit polynomial.

        Raises
        ------
        SlitFitError
            If ``simul_slit.wavelength`` is ``None`` or has a different shape
            from ``simul_slit.data``.
        SlitFitError
            If the number of valid pixels (after masking) is smaller than
            ``self.max_degree + 1``.

        Notes
        -----
        A pixel is included in the fit only when all of the following hold:

        * ``simul_slit.data != 0``  (simulation has signal here)
        * ``simul_slit.wavelength > 0``  (wavelength was assigned)
        * ``np.isfinite(observed_slit.data)``  (observed pixel is valid)
        * ``np.isfinite(simul_slit.data)``  (simulated pixel is valid)
        * DQ bit 0 (DO_NOT_USE) is not set in ``observed_slit.dq``
        """
        data_masked, sim_masked, lam, lam_ref = _build_fit_arrays(observed_slit, simul_slit)

        n_valid = len(data_masked)
        if n_valid < self.max_degree + 1:
            raise SlitFitError(
                f"Only {n_valid} valid pixel(s) available for a degree-{self.max_degree} "
                f"polynomial fit (need at least {self.max_degree + 1}). Reduce the fitting degree "
                "or check the input data."
            )

        dlam = lam - lam_ref
        coeffs = np.zeros(self.max_degree + 1)
        current_poly_vals = np.zeros(n_valid)

        for k in range(self.max_degree + 1):
            residual = data_masked - sim_masked * current_poly_vals
            col = sim_masked * dlam**k
            col_sq = float(np.dot(col, col))
            if col_sq == 0:
                break
            # Because the model is linear in the coefficients, the optimal c_k can be found
            # by a simple projection of the residual onto the new column
            c_k = float(np.dot(col, residual)) / col_sq

            if k >= 1 and self.improvement_threshold is not None:
                # Compute residual sum of squares
                rss_before = float(np.dot(residual, residual))
                new_residual = residual - c_k * col
                rss_after = float(np.dot(new_residual, new_residual))
                relative_improvement = (
                    (rss_before - rss_after) / rss_before if rss_before > 0 else 0.0
                )
                if relative_improvement < self.improvement_threshold:
                    log.debug(
                        f"Degree-{k} fit rejected: relative RSS improvement "
                        f"{relative_improvement:.6f} < threshold {self.improvement_threshold:.6f}"
                    )
                    break

            coeffs[k] = c_k
            current_poly_vals = current_poly_vals + c_k * dlam**k

        return _make_poly(coeffs, lam_ref)


def apply_flam_to_slit(sim_data, wavelength, f_lam):
    """
    Apply a fitted spectral polynomial to simulated slit data.

    Evaluates ``f_lam`` on the wavelength grid, zeros pixels outside the
    simulation footprint (where ``sim_data == 0`` or ``wavelength <= 0``),
    and returns the scaled simulated data.

    Parameters
    ----------
    sim_data : array-like
        2-D simulated flux array.
    wavelength : array-like
        2-D wavelength array with the same shape as ``sim_data``.
    f_lam : callable
        Function returned by a SlitFitter object.

    Returns
    -------
    scaled : ndarray
        The rescaled data.
    """
    sim_data = np.asarray(sim_data)
    wavelength = np.asarray(wavelength)
    footprint = (sim_data != 0) & (wavelength > 0)
    poly_surface = np.where(footprint, f_lam(wavelength), 0.0)
    return sim_data * poly_surface
