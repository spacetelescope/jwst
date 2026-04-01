import numpy as np
import pytest
from numpy.testing import assert_allclose
from stdatamodels.jwst.datamodels import SlitModel

from jwst.wfss_contam.wavefit import fit_spectral_shape


def _make_slits(shape, wl_region, true_coeffs, noise_level=0.0, rng=None):
    """
    Return observed and simulated slit mocks with a known polynomial spectrum.

    Parameters
    ----------
    shape : tuple
        (nrows, ncols) of the slit arrays.
    wl_region : tuple
        (row_slice, col_slice, wl_min, wl_max) defining where wavelengths are
        assigned and the simulation has signal.
    true_coeffs : sequence
        Polynomial coefficients [c0, c1, ...] for p(λ - λ_ref).
    """
    rows, cols, wl_min, wl_max = wl_region
    wavelength = np.zeros(shape)
    n_wl_cols = cols.stop - cols.start
    wavelength[rows, cols] = np.linspace(wl_min, wl_max, n_wl_cols)

    sim_data = np.zeros(shape)
    sim_data[rows, cols] = 1.0  # flat spectrum

    lam_ref = float(np.median(wavelength[wavelength > 0]))
    dlam = np.where(wavelength > 0, wavelength - lam_ref, 0.0)
    poly = sum(c * dlam**k for k, c in enumerate(true_coeffs))
    poly = np.where(sim_data != 0, poly, 0.0)
    obs_data = sim_data * poly

    if noise_level > 0.0:
        if rng is None:
            rng = np.random.default_rng(42)
        obs_data = obs_data + rng.normal(0, noise_level, shape)

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(shape, dtype=np.uint32)

    simul_slit = SlitModel()
    simul_slit.data = sim_data
    simul_slit.wavelength = wavelength

    return observed_slit, simul_slit, lam_ref


@pytest.mark.parametrize(
    "degree,true_coeffs",
    [
        (0, [2.5]),
        (1, [1.0, 0.5]),
        (2, [1.0, 0.5, -0.3]),
    ],
)
def test_fit_recovers_true_shape(degree, true_coeffs):
    """fit_spectral_shape scaled simulation reproduces observed data to within the noise level."""
    noise_level = 0.01
    shape = (20, 50)
    wl_region = (slice(5, 15), slice(5, 45), 1.7, 2.3)
    observed_slit, simul_slit, lam_ref = _make_slits(
        shape, wl_region, true_coeffs, noise_level=noise_level
    )

    scaled_sim, _, _ = fit_spectral_shape(observed_slit, simul_slit, degree=degree)
    assert_allclose(scaled_sim, observed_slit.data, atol=5 * noise_level)


def test_dq_mask_excludes_bad_pixels():
    """Pixels with DO_NOT_USE set (bit 0) are excluded from the fit; scaled simulation
    reproduces the good observed pixels to within the noise level."""
    noise_level = 0.01
    shape = (20, 50)
    wl_region = (slice(5, 15), slice(5, 45), 1.7, 2.3)
    # True spectrum is constant = 2.0
    observed_slit, simul_slit, _ = _make_slits(shape, wl_region, [2.0], noise_level=noise_level)

    # Corrupt some pixels in observed but mark them as bad
    corrupted = observed_slit.data.copy()
    bad_rows = slice(5, 7)
    corrupted[bad_rows, 5:45] = 999.0
    observed_slit.data = corrupted

    dq = np.zeros(shape, dtype=np.uint32)
    dq[bad_rows, 5:45] = 1  # DO_NOT_USE
    observed_slit.dq = dq

    scaled_sim, _, _ = fit_spectral_shape(observed_slit, simul_slit, degree=0)

    inside = simul_slit.data != 0
    good = (observed_slit.dq & 1) == 0
    assert_allclose(
        scaled_sim[inside & good], observed_slit.data[inside & good], atol=5 * noise_level
    )


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
        fit_spectral_shape(observed_slit, simul_slit, degree=2)


def test_raises_on_missing_wavelength():
    """Raises ValueError when simul_slit.wavelength is None."""
    sim_data = np.ones((5, 5))
    obs_data = np.ones((5, 5))

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros((5, 5), dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = sim_data
    simul_slit.wavelength = None

    with pytest.raises(ValueError, match="wavelength"):
        fit_spectral_shape(observed_slit, simul_slit)


def test_raises_on_wavelength_shape_mismatch():
    """Raises ValueError when wavelength shape differs from simul_slit.data shape."""
    sim_data = np.ones((5, 5))
    obs_data = np.ones((5, 5))

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros((5, 5), dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = sim_data
    simul_slit.wavelength = np.ones((3, 3))  # wrong shape

    with pytest.raises(ValueError, match="wavelength"):
        fit_spectral_shape(observed_slit, simul_slit)
