"""
Iterative fitting of a polynomial wavelength dependence to the data.

Idea is to take the simulated slit of each source, which initially assumes
a flat spectrum, and add in a polynomial wavelength dependence using slit.wavelength.
Then diff the simulation with the data for that slit and find the residual.
Finally optimize the polynomial coefficients to minimize the residuals.
"""

import numpy as np


def fit_spectral_shape(observed_slit, simul_slit, degree=2):
    """
    Fit a polynomial spectral shape to a simulated slit.

    The initially flat-spectrum ``simul_slit.data`` is scaled pixel-by-pixel
    by a polynomial in wavelength::

        scaled[y, x] = simul[y, x] * p(λ[y, x])
        p(λ) = Σ_k  c_k * (λ - λ_ref)^k

    The coefficients are chosen via linear least squares to minimize the
    squared residuals between the scaled simulation and the observed data::

        min  ||observed.data  -  simul.data * p(simul.wavelength)||²

    Because ``p`` is linear in the coefficients, the design matrix columns
    are simply ``simul.data * (λ - λ_ref)^k``, so no iterative solver is
    needed.

    Parameters
    ----------
    observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        The calibrated observed 2-D spectral cutout.
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        The flat-spectrum simulated cutout for the same source/order.
    degree : int, optional
        Degree of the fitting polynomial.

    Returns
    -------
    scaled_sim : ndarray
        The simulated data with the wavelength solution applied, same shape as ``simul_slit.data``.
        Zero where ``simul.data == 0`` or ``wavelength == 0``.

    Notes
    -----
    Pixels are included in the fit only when **all** of the following hold:

    * ``simul_slit.data != 0``  (simulation has signal here)
    * ``simul_slit.wavelength > 0``  (wavelength was assigned)
    * ``np.isfinite(observed_slit.data)``  (observed pixel is valid)
    * ``observed_slit.dq`` is not set to DO_NOT_USE

    If too few valid pixels remain for the requested polynomial degree, a
    ``ValueError`` is raised.
    """
    data = observed_slit.data
    sim = simul_slit.data
    wavelength = simul_slit.wavelength

    if wavelength is None or wavelength.shape != sim.shape:
        raise ValueError(
            "simul_slit.wavelength must be a 2-D array with the same shape as simul_slit.data"
        )

    # Build pixel validity mask
    dq = np.asarray(observed_slit.dq, dtype=np.uint32)
    mask = (
        (sim != 0)  # simulation has signal here
        & (wavelength > 0)  # wavelength was assigned
        & np.isfinite(data)  # observed pixel is valid
        & np.isfinite(sim)  # simulated pixel is valid
        & ((dq & 1) == 0)  # dq bit 0 is DO_NOT_USE
    )

    n_valid = int(np.sum(mask))
    if n_valid < degree + 1:
        raise ValueError(
            f"Only {n_valid} valid pixel(s) available for a degree-{degree} polynomial fit "
            f"(need at least {degree + 1}). Reduce the fitting degree or check the input data."
        )

    data_masked = data[mask]
    sim_masked = sim[mask]
    lam = wavelength[mask]

    # Center wavelengths around zero for numerical stability
    lam_ref = float(np.median(lam))
    dlam = lam - lam_ref

    # Design matrix: column k is  s * (λ - λ_ref)^k
    # coeffs in result are: p(λ) = sum(coeffs[k] * (λ - lam_ref)**k  for k in range(degree+1))
    design_matrix = np.column_stack([sim_masked * dlam**k for k in range(degree + 1)])
    coeffs, _residuals, _rank, _sv = np.linalg.lstsq(design_matrix, data_masked, rcond=None)

    # Evaluate the polynomial over the full 2-D footprint of the simulation.
    # Pixels outside the simulation footprint stay at zero.
    sim_footprint = (sim != 0) & (wavelength > 0)
    dlam_2d = np.where(sim_footprint, wavelength - lam_ref, 0.0)
    poly_surface = np.zeros(sim.shape, dtype=float)
    for k in range(degree + 1):
        poly_surface += coeffs[k] * dlam_2d**k
    poly_surface = np.where(sim_footprint, poly_surface, 0.0)

    scaled_sim = sim * poly_surface
    return scaled_sim
