"""Fit a spectral shape to a WFSS slit."""

import numpy as np

__all__ = ["SlitPolynomialFitter", "apply_flam_to_slit"]


class SlitFitError(Exception):
    """Base class for exceptions in this module."""

    pass


class SlitPolynomialFitter:
    """
    Fit a polynomial spectral shape to observed and simulated slits.

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

        The coefficients are chosen via linear least squares to minimize the
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
            A function ``flam(wavelength)`` that evaluates the best-fit
            polynomial $p(λ)$ at every point in a wavelength array and
            returns an array of the same shape.

        Raises
        ------
        ValueError
            If ``simul_slit.wavelength`` is ``None`` or has a different shape
            from ``simul_slit.data``.
        ValueError
            If the number of valid pixels (after masking) is smaller than
            ``self.degree + 1``.

        Notes
        -----
        A pixel is included in the fit only when **all** of the following hold:

        * ``simul_slit.data != 0``  (simulation has signal here)
        * ``simul_slit.wavelength > 0``  (wavelength was assigned)
        * ``np.isfinite(observed_slit.data)``  (observed pixel is valid)
        * ``np.isfinite(simul_slit.data)``  (simulated pixel is valid)
        * DQ bit 0 (DO_NOT_USE) is not set in ``observed_slit.dq``

        Before fitting, a single 3-sigma clip is applied to remove outliers in
        the observed data residuals from the initial flat-spectrum ratio.
        """
        degree = self.degree
        data = np.asarray(observed_slit.data)
        sim = np.asarray(simul_slit.data)
        wavelength = np.asarray(simul_slit.wavelength)

        if wavelength is None or wavelength.shape != sim.shape:
            raise SlitFitError(
                "simul_slit.wavelength must be a 2-D array with the same shape as simul_slit.data"
            )

        # Build pixel validity mask
        dq = np.asarray(observed_slit.dq, dtype=np.uint32)
        base_mask = (
            (sim != 0)  # simulation has signal here
            & (wavelength > 0)  # wavelength was assigned
            & np.isfinite(data)  # observed pixel is valid
            & np.isfinite(sim)  # simulated pixel is valid
            & ((dq & 1) == 0)  # dq bit 0 is DO_NOT_USE
        )

        n_valid = int(np.sum(base_mask))
        if n_valid < degree + 1:
            raise SlitFitError(
                f"Only {n_valid} valid pixel(s) available for a degree-{degree} polynomial fit "
                f"(need at least {degree + 1}). Reduce the fitting degree or check the input data."
            )

        data_masked = data[base_mask]
        sim_masked = sim[base_mask]
        lam = wavelength[base_mask]

        # Center wavelengths around zero for numerical stability
        lam_ref = float(np.median(lam))
        dlam = lam - lam_ref

        # Design matrix: column k is  s * (λ - λ_ref)^k
        design_matrix = np.column_stack([sim_masked * dlam**k for k in range(degree + 1)])
        coeffs, *_ = np.linalg.lstsq(design_matrix, data_masked, rcond=None)

        def f_lam(wavelength):
            # Return should be an f(wavelength) so we can apply it to multiple slits
            wavelength = np.asarray(wavelength)
            dlam = wavelength - lam_ref
            poly_vals = np.zeros(wavelength.shape, dtype=float)
            for k in range(degree + 1):
                poly_vals += coeffs[k] * dlam**k
            return poly_vals

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
        Function returned by a SlitFitter

    Returns
    -------
    scaled : ndarray
        ``sim_data`` scaled by the polynomial, zeroed outside the footprint.
    """
    sim_data = np.asarray(sim_data)
    wavelength = np.asarray(wavelength)
    footprint = (sim_data != 0) & (wavelength > 0)
    poly_surface = np.where(footprint, f_lam(wavelength), 0.0)
    return sim_data * poly_surface
