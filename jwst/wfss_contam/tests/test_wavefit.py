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
    """
    Make a pair of simulated and observed SlitModel objects.

    Observed data has a sinusoidal spectral shape plus noise; simulated data is initially flat.

    Returns
    -------
    observed_slit
        Mock observed slit with a sinusoidal spectral shape.
    simul_slit
        Simulated slit with a flat spectrum.
    wavelength
        Wavelengths used as input to create the mock observed wavelength-dependent data.
    """
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
    # odd nor purely even about lam_ref.
    # Otherwise some basis functions might not contribute to the fit.
    period = wl_max - wl_min
    sinusoid = np.sin(np.pi * wavelength / period)
    noise_level = 0.1
    rng = np.random.default_rng(42)
    noise = rng.normal(0, noise_level, shape)
    obs_data = np.where(sim_data != 0, sim_data + sinusoid + noise, 0.0)

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(shape, dtype=np.uint32)

    simul_slit = SlitModel()
    simul_slit.data = sim_data
    return observed_slit, simul_slit, wavelength


def test_more_basis_images_reduces_residuals(make_slits):
    """
    Adding more ``fluxmodel_N`` basis images progressively reduces fit residuals.

    The observed data follow a sinusoidal spectral shape.  The basis images are
    sim * (λ - λ_ref)^k (k = 0, 1, 2, ...), so each extra term adds one more
    Taylor-series contribution that helps approximate the sinusoid.  The RMS
    between the fitted simulation and the observed data should decrease as
    more basis images are used.
    """
    observed_slit, simul_slit, wavelength = make_slits
    obs_data = observed_slit.data
    sim_data = simul_slit.data
    inside = sim_data != 0

    lam_ref = float(np.median(wavelength[inside]))
    dlam = np.where(inside, wavelength - lam_ref, 0.0)

    max_n_terms = 5
    rms_values = []
    for n_terms in range(1, max_n_terms + 1):
        test_simul = SlitModel()
        test_simul.data = sim_data.copy()
        for k in range(1, n_terms):
            setattr(test_simul, f"fluxmodel_{k}", sim_data * dlam**k)

        coeffs = fit_slit_by_basis_images(observed_slit, test_simul)
        fitted = apply_basis_coeffs(test_simul, coeffs)
        rms = np.sqrt(np.mean((fitted[inside] - obs_data[inside]) ** 2))
        rms_values.append(rms)

    for i in range(len(rms_values) - 1):
        assert rms_values[i + 1] < rms_values[i]


def _make_basis_slit(sim_data, wavelength, max_order=1):
    """Build a simul_slit with max_order fluxmodel_k attributes."""
    inside = sim_data != 0
    lam_ref = float(np.median(wavelength[inside]))
    dlam = np.where(inside, wavelength - lam_ref, 0.0)
    slit = SlitModel()
    slit.data = sim_data.copy()
    for k in range(1, max_order + 1):
        setattr(slit, f"fluxmodel_{k}", sim_data * dlam**k)
    return slit


def test_fit_recovers_exact_coefficients(make_slits):
    """Test that fit_slit_by_basis_images recovers known coefficients exactly for noise-free case."""
    _, simul_slit, wavelength = make_slits
    sim_data = simul_slit.data

    c0_true, c1_true = 1.0, 1.5
    simul = _make_basis_slit(sim_data, wavelength, max_order=1)

    # Build observed as an exact linear combination of the two basis images
    obs_data = c0_true * simul.data + c1_true * simul.fluxmodel_1
    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(obs_data.shape, dtype=np.uint32)

    coeffs = fit_slit_by_basis_images(observed_slit, simul)
    assert coeffs is not None
    assert_allclose(coeffs, [c0_true, c1_true], rtol=1e-6)


def test_fit_flat_only_no_fluxmodel(make_slits):
    """With no fluxmodel_N attributes, fit returns a single scalar coefficient."""
    _, simul_slit, _ = make_slits
    sim_data = simul_slit.data

    scale = 1.05
    obs_data = scale * sim_data
    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(obs_data.shape, dtype=np.uint32)

    simul = SlitModel()
    simul.data = sim_data.copy()

    coeffs = fit_slit_by_basis_images(observed_slit, simul)
    assert coeffs is not None
    assert coeffs.shape == (1,)
    assert_allclose(coeffs[0], scale, rtol=1e-6)


def test_fit_dq_mask_excludes_bad_pixels(make_slits):
    """Test that masking of bad pixels is happening."""
    _, simul_slit, wavelength = make_slits
    sim_data = simul_slit.data

    c0_true, c1_true = 1.05, -0.5
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
    assert coeffs is not None
    assert_allclose(coeffs, [c0_true, c1_true], rtol=1e-6)


def test_large_errors_downweight_bad_pixels(make_slits):
    """Test that large error values on bad pixels downweight them and we can recover coefficients."""
    _, simul_slit, wavelength = make_slits
    sim_data = simul_slit.data

    c0_true, c1_true = 1.05, -0.5
    simul = _make_basis_slit(sim_data, wavelength, max_order=1)
    clean_obs = c0_true * simul.data + c1_true * simul.fluxmodel_1

    # Corrupt a patch, leave unmasked, set errors to be very large
    corrupted = clean_obs.copy()
    corrupted[5:8, 10:40] = 999.0
    err = np.ones(corrupted.shape)
    err[5:8, 10:40] = 1e10

    observed = SlitModel()
    observed.data = corrupted
    observed.dq = observed.get_default("dq")
    observed.err = err

    coeffs = fit_slit_by_basis_images(observed, simul)
    assert coeffs is not None
    assert_allclose(coeffs, [c0_true, c1_true], rtol=1e-6)


def test_fit_raises_when_all_errors_nonfinite(make_slits):
    """Test SlitFitError raised when the err array contains no finite positive values."""
    _, simul_slit, wavelength = make_slits
    sim_data = simul_slit.data

    simul = _make_basis_slit(sim_data, wavelength, max_order=1)
    obs_data = 1.0 * simul.data + 0.5 * simul.fluxmodel_1

    observed = SlitModel()
    observed.data = obs_data
    observed.dq = np.zeros(obs_data.shape, dtype=np.uint32)
    observed.err = np.full(obs_data.shape, np.nan)

    with pytest.raises(SlitFitError, match="finite positive error"):
        fit_slit_by_basis_images(observed, simul)


def test_fit_raises_too_few_valid_pixels():
    """Test SlitFitError raised when valid pixels < number of basis terms."""
    shape = (5, 5)
    sim_data = np.zeros(shape)
    sim_data[2, 2] = 1.0  # only 1 valid pixel

    simul = SlitModel()
    simul.data = sim_data
    simul.fluxmodel_1 = sim_data * 0.5

    observed = SlitModel()
    observed.data = sim_data.copy()
    observed.dq = np.zeros(shape, dtype=np.uint32)

    with pytest.raises(SlitFitError, match="valid pixel"):
        fit_slit_by_basis_images(observed, simul)


def test_fit_returns_none_when_c0_far_from_unity(make_slits):
    """
    Test that fit_slit_by_basis_images returns None when c_0 is far from 1.

    For now this is a hard-coded 'bad fit' criterion.
    """
    _, simul_slit, wavelength = make_slits
    sim_data = simul_slit.data

    simul = _make_basis_slit(sim_data, wavelength, max_order=1)
    obs_data = 5.0 * simul.data
    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(obs_data.shape, dtype=np.uint32)

    result = fit_slit_by_basis_images(observed_slit, simul)
    assert result is None


def test_apply_basis_coeffs_correct_linear_combination(make_slits):
    """Test that apply_basis_coeffs computes linear combinations as expected."""
    _, simul_slit, wavelength = make_slits
    sim_data = simul_slit.data

    simul = _make_basis_slit(sim_data, wavelength, max_order=1)
    c0, c1 = 1.8, -0.3
    coeffs = np.array([c0, c1])

    result = apply_basis_coeffs(simul, coeffs)
    expected = c0 * np.asarray(simul.data) + c1 * np.asarray(simul.fluxmodel_1)
    assert_allclose(result, expected, rtol=1e-6)
