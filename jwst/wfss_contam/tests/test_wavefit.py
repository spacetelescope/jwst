import numpy as np
import pytest
from numpy.testing import assert_allclose
from stdatamodels.jwst.datamodels import SlitModel

from jwst.wfss_contam.wavefit import (
    SlitFitError,
    apply_basis_coeffs,
    fit_slit_by_basis_images,
)


@pytest.fixture
def make_slits():
    shape = (20, 100)
    rows = slice(5, 15)
    cols = slice(5, 95)
    wl_min, wl_max = 1.0, 3.0

    wavelength = np.zeros(shape)
    n_wl_cols = cols.stop - cols.start
    wavelength[rows, cols] = np.linspace(wl_min, wl_max, n_wl_cols)

    sim_data = np.zeros(shape)
    sim_data[rows, cols] = 1.0

    # Sinusoidal spectrum with a phase offset so the function is neither purely
    # odd nor purely even about lam_ref, ensuring every polynomial degree contributes.
    period = wl_max - wl_min  # one full cycle across the wavelength range
    sinusoid = np.sin(np.pi * wavelength / period)
    noise_level = 0.1
    rng = np.random.default_rng(42)
    noise = rng.normal(0, noise_level, shape)
    obs_data = np.where(sim_data != 0, sinusoid + noise, 0.0)

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(shape, dtype=np.uint32)

    simul_slit = SlitModel()
    simul_slit.data = sim_data
    simul_slit.wavelength = wavelength
    return observed_slit, simul_slit


def test_more_basis_images_reduces_residuals(make_slits):
    """
    Adding more fluxmodel_N basis images progressively reduces fit residuals.

    The observed data follow a sinusoidal spectral shape.  The basis images are
    sim * (λ - λ_ref)^k (k = 0, 1, 2, ...), so each extra term adds one more
    Taylor-series contribution that helps approximate the sinusoid.  The RMS
    between the fitted simulation and the observed data must strictly decrease
    as basis images are added.
    """
    observed_slit, simul_slit = make_slits
    obs_data = observed_slit.data
    sim_data = simul_slit.data
    wavelength = simul_slit.wavelength
    inside = sim_data != 0

    lam_ref = float(np.median(wavelength[inside]))
    dlam = np.where(inside, wavelength - lam_ref, 0.0)

    max_n_terms = 5
    rms_values = []
    for n_terms in range(1, max_n_terms + 1):
        test_simul = SlitModel()
        test_simul.data = sim_data.copy()
        test_simul.wavelength = wavelength.copy()
        for k in range(1, n_terms):
            setattr(test_simul, f"fluxmodel_{k}", sim_data * dlam**k)

        coeffs = fit_slit_by_basis_images(observed_slit, test_simul)
        fitted = apply_basis_coeffs(test_simul, coeffs)
        rms = np.sqrt(np.mean((fitted[inside] - obs_data[inside]) ** 2))
        rms_values.append(rms)

    for i in range(len(rms_values) - 1):
        assert rms_values[i + 1] < rms_values[i], (
            f"RMS did not decrease from {i + 1} basis term(s) ({rms_values[i]:.6f}) "
            f"to {i + 2} basis term(s) ({rms_values[i + 1]:.6f})"
        )


def _make_basis_slit(sim_data, wavelength, max_order=1):
    """Build a simul_slit with max_order fluxmodel_k attributes."""
    inside = sim_data != 0
    lam_ref = float(np.median(wavelength[inside]))
    dlam = np.where(inside, wavelength - lam_ref, 0.0)
    slit = SlitModel()
    slit.data = sim_data.copy()
    slit.wavelength = wavelength.copy()
    for k in range(1, max_order + 1):
        setattr(slit, f"fluxmodel_{k}", sim_data * dlam**k)
    return slit


def test_fit_recovers_exact_coefficients(make_slits):
    """fit_slit_by_basis_images recovers known coefficients exactly (no noise)."""
    _, simul_slit = make_slits
    sim_data = simul_slit.data
    wavelength = simul_slit.wavelength

    c0_true, c1_true = 3.0, 1.5
    simul = _make_basis_slit(sim_data, wavelength, max_order=1)

    # Build observed as an exact linear combination of the two basis images
    obs_data = c0_true * simul.data + c1_true * simul.fluxmodel_1
    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(obs_data.shape, dtype=np.uint32)

    coeffs = fit_slit_by_basis_images(observed_slit, simul)
    assert_allclose(coeffs, [c0_true, c1_true], rtol=1e-10)


def test_fit_flat_only_no_fluxmodel(make_slits):
    """With no fluxmodel_N attributes, fit returns a single scalar coefficient."""
    _, simul_slit = make_slits
    sim_data = simul_slit.data

    scale = 2.7
    obs_data = scale * sim_data
    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(obs_data.shape, dtype=np.uint32)

    simul = SlitModel()
    simul.data = sim_data.copy()
    simul.wavelength = simul_slit.wavelength.copy()

    coeffs = fit_slit_by_basis_images(observed_slit, simul)
    assert coeffs.shape == (1,)
    assert_allclose(coeffs[0], scale, rtol=1e-6)


def test_fit_dq_mask_excludes_bad_pixels(make_slits):
    """Coefficients are unaffected when bad pixels are masked."""
    _, simul_slit = make_slits
    sim_data = simul_slit.data
    wavelength = simul_slit.wavelength

    c0_true, c1_true = 2.0, -0.5
    simul = _make_basis_slit(sim_data, wavelength, max_order=1)
    clean_obs = c0_true * simul.data + c1_true * simul.fluxmodel_1

    # Corrupt a patch and mask it
    corrupted = clean_obs.copy()
    corrupted[5:8, 10:40] = 999.0
    dq = np.zeros(corrupted.shape, dtype=np.uint32)
    dq[5:8, 10:40] = 1

    observed_masked = SlitModel()
    observed_masked.data = corrupted
    observed_masked.dq = dq

    coeffs = fit_slit_by_basis_images(observed_masked, simul)
    assert_allclose(coeffs, [c0_true, c1_true], rtol=1e-10)


def test_fit_raises_too_few_valid_pixels():
    """SlitFitError is raised when valid pixels < number of basis terms."""
    shape = (5, 5)
    sim_data = np.zeros(shape)
    sim_data[2, 2] = 1.0  # only 1 valid pixel

    simul = SlitModel()
    simul.data = sim_data
    simul.wavelength = np.zeros(shape)
    simul.wavelength[2, 2] = 2.0
    simul.fluxmodel_1 = sim_data * 0.5

    observed = SlitModel()
    observed.data = sim_data.copy()
    observed.dq = np.zeros(shape, dtype=np.uint32)

    with pytest.raises(SlitFitError, match="valid pixel"):
        fit_slit_by_basis_images(observed, simul)


def test_apply_basis_coeffs_correct_linear_combination(make_slits):
    """apply_basis_coeffs returns c0*data + c1*fluxmodel_1 exactly."""
    _, simul_slit = make_slits
    sim_data = simul_slit.data
    wavelength = simul_slit.wavelength

    simul = _make_basis_slit(sim_data, wavelength, max_order=1)
    c0, c1 = 1.8, -0.3
    coeffs = np.array([c0, c1])

    result = apply_basis_coeffs(simul, coeffs)
    expected = c0 * np.asarray(simul.data, dtype=float) + c1 * np.asarray(
        simul.fluxmodel_1, dtype=float
    )
    assert_allclose(result, expected, rtol=1e-6)


def test_apply_basis_coeffs_raises_on_wrong_length(make_slits):
    """apply_basis_coeffs raises ValueError when len(coeffs) != number of basis images."""
    _, simul_slit = make_slits
    simul = _make_basis_slit(simul_slit.data, simul_slit.wavelength, max_order=1)
    # 2 basis images but 3 coefficients
    with pytest.raises(ValueError):
        apply_basis_coeffs(simul, np.array([1.0, 2.0, 3.0]))
