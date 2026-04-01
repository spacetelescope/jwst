import numpy as np
import pytest
from numpy.testing import assert_allclose
from stdatamodels.jwst.datamodels import SlitModel

from jwst.wfss_contam.wavefit import SlitPolynomialFitter


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
    noise_level = 0.01
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


def test_higher_degree_reduces_residuals_for_sinusoidal_spectrum(make_slits):
    """Higher-degree polynomial fits produce smaller residuals on a sinusoidal spectrum."""
    observed_slit, simul_slit = make_slits
    obs_data = observed_slit.data
    sim_data = simul_slit.data

    inside = sim_data != 0
    rms_values = []
    for degree in range(5):
        scaled_sim = SlitPolynomialFitter(degree=degree)(observed_slit, simul_slit)
        rms = np.sqrt(np.mean((scaled_sim[inside] - obs_data[inside]) ** 2))
        rms_values.append(rms)

    # Each additional degree should strictly improve the fit
    for i in range(len(rms_values) - 1):
        assert rms_values[i + 1] < rms_values[i], (
            f"RMS did not decrease from degree {i} ({rms_values[i]:.6f}) "
            f"to degree {i + 1} ({rms_values[i + 1]:.6f})"
        )


def test_dq_mask_excludes_bad_pixels(make_slits):
    """Ensure pixels with DO_NOT_USE set are excluded from the fit."""
    observed_slit, simul_slit = make_slits

    # Corrupt some pixels in observed but give them DO_NOT_USE flag
    corrupted = observed_slit.data.copy()
    bad_rows = slice(5, 7)
    corrupted[bad_rows, 5:45] = 999.0
    observed_slit.data = corrupted

    dq = np.zeros(observed_slit.data.shape, dtype=np.uint32)
    dq[bad_rows, 5:45] = 1
    observed_slit.dq = dq

    # fit with dq flags: should ignore bad pixels
    fitter = SlitPolynomialFitter(degree=5)
    scaled_sim_masked = fitter(observed_slit, simul_slit)

    inside = simul_slit.data != 0
    good = (observed_slit.dq & 1) == 0
    assert_allclose(scaled_sim_masked[inside & good], observed_slit.data[inside & good], atol=0.05)

    # fit without flags: should give a much worse fit
    observed_slit_unmasked = observed_slit.copy()
    observed_slit_unmasked.dq = np.zeros_like(dq)
    scaled_sim_unmasked = fitter(observed_slit_unmasked, simul_slit)

    rms_masked = np.sqrt(
        np.mean((scaled_sim_masked[inside & good] - observed_slit.data[inside & good]) ** 2)
    )
    rms_unmasked = np.sqrt(
        np.mean((scaled_sim_unmasked[inside & good] - observed_slit.data[inside & good]) ** 2)
    )
    assert rms_unmasked > 100 * rms_masked


def test_raises_when_too_few_pixels_for_degree():
    """Raises ValueError when valid pixel count < degree+1."""
    shape = (5, 5)
    # Only one valid pixel
    sim_data = np.zeros(shape)
    sim_data[2, 2] = 1.0
    wavelength = np.zeros(shape)
    wavelength[2, 2] = 2.0
    obs_data = np.zeros(shape)
    obs_data[2, 2] = 3.0

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(shape, dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = sim_data
    simul_slit.wavelength = wavelength

    # degree=2 requires 3 pixels; only 1 available
    with pytest.raises(ValueError, match="valid pixel"):
        SlitPolynomialFitter()(observed_slit, simul_slit)


def test_error_missing_wavelength():
    """Raises ValueError when simul_slit.wavelength is None."""
    observed_slit = SlitModel()
    observed_slit.data = np.ones((5, 5))
    observed_slit.dq = np.zeros((5, 5), dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = np.ones((5, 5))
    simul_slit.wavelength = None

    with pytest.raises(ValueError, match="wavelength"):
        SlitPolynomialFitter()(observed_slit, simul_slit)


def test_error_wavelength_shape_mismatch():
    """Raises ValueError when wavelength shape differs from simul_slit.data shape."""
    observed_slit = SlitModel()
    observed_slit.data = np.ones((5, 5))
    observed_slit.dq = np.zeros((5, 5), dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = np.ones((5, 5))
    simul_slit.wavelength = np.ones((3, 3))  # wrong shape

    with pytest.raises(ValueError, match="wavelength"):
        SlitPolynomialFitter()(observed_slit, simul_slit)
