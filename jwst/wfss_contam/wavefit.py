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
    Iteratively fit polynomials of increasing degree to a WFSS slit.

    Fits are performed successively: after each step the current simulated
    spectrum is multiplied by the fitted polynomial before the next fit is
    performed.  For example, ``degree_sequence=[1, 5]`` first fits a
    degree-1 polynomial to capture the broad spectral shape, then fits a
    degree-5 polynomial to the residuals, and returns their product as the
    final spectral shape function.

    Parameters
    ----------
    degree_sequence : list of int
        Sequence of polynomial degrees to fit in order.
    improvement_threshold : float or None, optional
        Minimum required relative improvement in the residual sum-of-squares
        (RSS) for each successive fit step to be accepted.  Computed as
        ``(RSS_before - RSS_after) / RSS_before``.  If a fit step improves
        the RSS by less than this fraction, that polynomial is rejected and
        iteration stops early.  ``None`` (default) disables the check and
        always applies every step in ``degree_sequence``.
    """

    def __init__(self, degree_sequence, improvement_threshold=None):
        self.degree_sequence = list(degree_sequence)
        self.improvement_threshold = improvement_threshold

    def __call__(self, observed_slit, simul_slit):
        """
        Scale the simulation by a product of successively fitted polynomials.

        At each step the current simulated spectrum (initialised to
        ``simul_slit.data``) is multiplied by the polynomial found at the previous step,
        and a new polynomial ``p_i`` is fitted to::

            min  ||observed.data  -  current_sim * p_i(λ)||²

        The function returned is the product of all per-step polynomials::

            f(λ) = p_1(λ) · p_2(λ) · … · p_n(λ)

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
            A function ``f_lam(wavelength)`` that evaluates the combined
            best-fit polynomial spectral shape.

        Raises
        ------
        SlitFitError
            If ``simul_slit.wavelength`` is ``None`` or has a different shape
            from ``simul_slit.data``.
        SlitFitError
            If the number of valid pixels (after masking) is smaller than
            ``max(self.degree_sequence) + 1``.

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
        max_degree = max(self.degree_sequence)
        if n_valid < max_degree + 1:
            raise SlitFitError(
                f"Only {n_valid} valid pixel(s) available for a degree-{max_degree} polynomial fit "
                f"(need at least {max_degree + 1}). Reduce the fitting degree "
                "or check the input data."
            )

        dlam = lam - lam_ref
        f_lam_factors = []
        current_sim_masked = sim_masked.copy()
        for degree in self.degree_sequence:
            design_matrix = np.column_stack(
                [current_sim_masked * dlam**k for k in range(degree + 1)]
            )
            coeffs, *_ = np.linalg.lstsq(design_matrix, data_masked, rcond=None)
            f_k = _make_poly(coeffs.copy(), lam_ref)

            if self.improvement_threshold is not None:
                rss_before = float(np.sum((data_masked - current_sim_masked) ** 2))
                rss_after = float(np.sum((data_masked - current_sim_masked * f_k(lam)) ** 2))
                relative_improvement = (
                    (rss_before - rss_after) / rss_before if rss_before > 0 else 0.0
                )
                if relative_improvement < self.improvement_threshold:
                    log.debug(
                        f"Degree-{degree} fit rejected: relative RSS improvement "
                        f"{relative_improvement:.3f} < threshold {self.improvement_threshold:.3f}"
                    )
                    break

            f_lam_factors.append(f_k)
            current_sim_masked = current_sim_masked * f_k(lam)

        def f_lam(wavelength, _factors=f_lam_factors):
            wavelength = np.asarray(wavelength)
            result = np.ones(wavelength.shape, dtype=float)
            for f in _factors:
                result *= f(wavelength)
            return result

        return f_lam


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
