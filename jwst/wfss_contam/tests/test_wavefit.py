import numpy as np
import pytest
from numpy.testing import assert_allclose
from stdatamodels.jwst.datamodels import SlitModel

from jwst.wfss_contam.wavefit import (
    SlitFitError,
    SlitIterativePolynomialFitter,
    SlitPolynomialFitter,
    _build_fit_arrays,
    apply_flam_to_slit,
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


def test_higher_degree_reduces_residuals_for_sinusoidal_spectrum(make_slits):
    """Higher-degree polynomial fits produce smaller residuals on a sinusoidal spectrum."""
    observed_slit, simul_slit = make_slits
    obs_data = observed_slit.data
    sim_data = simul_slit.data

    inside = sim_data != 0
    rms_values = []
    for degree in range(5):
        poly_fn = SlitPolynomialFitter(degree=degree)(observed_slit, simul_slit)
        scaled_sim = apply_flam_to_slit(simul_slit.data, simul_slit.wavelength, poly_fn)
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
    poly_fn_masked = fitter(observed_slit, simul_slit)
    scaled_sim_masked = apply_flam_to_slit(simul_slit.data, simul_slit.wavelength, poly_fn_masked)

    inside = simul_slit.data != 0
    good = (observed_slit.dq & 1) == 0
    # atol is 5x RMS
    assert_allclose(scaled_sim_masked[inside & good], observed_slit.data[inside & good], atol=0.5)

    # fit without flags: should give a much worse fit
    observed_slit_unmasked = observed_slit.copy()
    observed_slit_unmasked.dq = np.zeros_like(dq)
    poly_fn_unmasked = fitter(observed_slit_unmasked, simul_slit)
    scaled_sim_unmasked = apply_flam_to_slit(
        simul_slit.data, simul_slit.wavelength, poly_fn_unmasked
    )

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
    with pytest.raises(SlitFitError, match="valid pixel"):
        SlitPolynomialFitter()(observed_slit, simul_slit)


def test_error_missing_wavelength():
    """Raises ValueError when simul_slit.wavelength is None."""
    observed_slit = SlitModel()
    observed_slit.data = np.ones((5, 5))
    observed_slit.dq = np.zeros((5, 5), dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = np.ones((5, 5))
    simul_slit.wavelength = None

    with pytest.raises(SlitFitError, match="wavelength"):
        SlitPolynomialFitter()(observed_slit, simul_slit)


def test_error_wavelength_shape_mismatch():
    """Raises ValueError when wavelength shape differs from simul_slit.data shape."""
    observed_slit = SlitModel()
    observed_slit.data = np.ones((5, 5))
    observed_slit.dq = np.zeros((5, 5), dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = np.ones((5, 5))
    simul_slit.wavelength = np.ones((3, 3))  # wrong shape

    with pytest.raises(SlitFitError, match="wavelength"):
        SlitPolynomialFitter()(observed_slit, simul_slit)


def test_iterative_single_element(make_slits):
    """SlitIterativePolynomialFitter([d]) produces the same f_lam as SlitPolynomialFitter(d)."""
    observed_slit, simul_slit = make_slits
    wl = simul_slit.wavelength

    f_single = SlitPolynomialFitter(degree=3)(observed_slit, simul_slit)
    f_iter = SlitIterativePolynomialFitter([3])(observed_slit, simul_slit)

    assert_allclose(f_single(wl), f_iter(wl), rtol=1e-12)


def test_iterative_fit_reduces_residuals_relative_to_first_degree(make_slits):
    """Iterative [low, high] fit gives at least as small residuals as low degree alone."""
    observed_slit, simul_slit = make_slits
    obs_data = observed_slit.data
    sim_data = simul_slit.data
    inside = sim_data != 0

    f_low = SlitPolynomialFitter(degree=1)(observed_slit, simul_slit)
    f_iter = SlitIterativePolynomialFitter([1, 5])(observed_slit, simul_slit)

    def rms(f):
        scaled = apply_flam_to_slit(sim_data, simul_slit.wavelength, f)
        return np.sqrt(np.mean((scaled[inside] - obs_data[inside]) ** 2))

    assert rms(f_iter) < rms(f_low)


def test_iterative_fit_function_product(make_slits):
    """
    The composed f_lam from SlitIterativePolynomialFitter([d1, d2]) equals the product
    of the two individual polynomials evaluated on the same wavelength grid.
    """
    observed_slit, simul_slit = make_slits
    wl = simul_slit.wavelength

    fitter = SlitIterativePolynomialFitter([1, 3])
    f_composed = fitter(observed_slit, simul_slit)

    data_m, sim_m, lam, lam_ref = _build_fit_arrays(observed_slit, simul_slit)
    dlam = lam - lam_ref

    # --- step 1: degree-1 ---
    dm1 = np.column_stack([sim_m * dlam**k for k in range(2)])
    c1, *_ = np.linalg.lstsq(dm1, data_m, rcond=None)
    p1 = sum(c1[k] * (wl - lam_ref) ** k for k in range(2))
    sim_m2 = sim_m * sum(c1[k] * dlam**k for k in range(2))

    # --- step 2: degree-3 ---
    dm2 = np.column_stack([sim_m2 * dlam**k for k in range(4)])
    c2, *_ = np.linalg.lstsq(dm2, data_m, rcond=None)
    p2 = sum(c2[k] * (wl - lam_ref) ** k for k in range(4))

    expected = p1 * p2
    assert_allclose(f_composed(wl), expected, rtol=1e-5)


def test_raises_when_too_few_pixels_for_degree_iterative():
    """Raises SlitFitError when valid pixels < max(degree_sequence) + 1."""
    shape = (5, 5)
    # 3 valid pixels — enough for degree 2, but not for degree 3
    sim_data = np.zeros(shape)
    wavelength = np.zeros(shape)
    obs_data = np.zeros(shape)
    for i, (r, c) in enumerate([(2, 1), (2, 2), (2, 3)]):
        sim_data[r, c] = 1.0
        wavelength[r, c] = 1.0 + i * 0.5
        obs_data[r, c] = 2.0

    observed_slit = SlitModel()
    observed_slit.data = obs_data
    observed_slit.dq = np.zeros(shape, dtype=np.uint32)
    simul_slit = SlitModel()
    simul_slit.data = sim_data
    simul_slit.wavelength = wavelength

    with pytest.raises(SlitFitError, match="valid pixel"):
        SlitIterativePolynomialFitter([1, 3])(observed_slit, simul_slit)


def test_improvement_threshold_zero_matches_no_threshold(make_slits):
    """improvement_threshold=0.0 accepts every step, producing the same result as None."""
    observed_slit, simul_slit = make_slits
    wl = simul_slit.wavelength

    f_none = SlitIterativePolynomialFitter([1, 3, 5])(observed_slit, simul_slit)
    f_zero = SlitIterativePolynomialFitter([1, 3, 5], improvement_threshold=0.0)(
        observed_slit, simul_slit
    )

    assert_allclose(f_zero(wl), f_none(wl), rtol=1e-12)


def test_improvement_threshold_blocks_all_steps(make_slits):
    """An impossibly high threshold blocks every step; f_lam returns ones everywhere."""
    observed_slit, simul_slit = make_slits
    wl = simul_slit.wavelength
    inside = simul_slit.data != 0

    f_blocked = SlitIterativePolynomialFitter([1, 3], improvement_threshold=2.0)(
        observed_slit, simul_slit
    )
    f_full = SlitIterativePolynomialFitter([1, 3])(observed_slit, simul_slit)

    assert_allclose(f_blocked(wl), np.ones_like(wl))


def test_improvement_threshold_basic(make_slits):
    """
    A moderate threshold stops before a very high-degree step but still fits well.

    Degree sequence [1, 5, 50] is used.  The first two steps capture the sinusoidal
    shape; by the time the degree-50 step is attempted, residuals are already near the
    noise floor so relative RSS improvement falls below the threshold and iteration
    stops.  The result should be:
      - better than no fit (rms < rms of flat simulation vs observation), and
      - worse than the unconstrained full sequence (rms > rms of full fit).
    """
    observed_slit, simul_slit = make_slits
    inside = simul_slit.data != 0

    def rms(f):
        scaled = apply_flam_to_slit(simul_slit.data, simul_slit.wavelength, f)
        return np.sqrt(np.mean((scaled[inside] - observed_slit.data[inside]) ** 2))

    # Baseline: flat simulation with no spectral correction
    rms_no_fit = np.sqrt(np.mean((simul_slit.data[inside] - observed_slit.data[inside]) ** 2))

    degrees = np.arange(1, 8)
    f_full = SlitIterativePolynomialFitter(degrees)(observed_slit, simul_slit)
    f_stopped = SlitIterativePolynomialFitter(degrees, improvement_threshold=0.3)(
        observed_slit, simul_slit
    )

    rms_full = rms(f_full)
    rms_stopped = rms(f_stopped)

    # Early-stopped fit is still a meaningful improvement over no correction
    assert rms_stopped < rms_no_fit
    # But the unconstrained fit (which runs degree-7) is better
    assert rms_stopped > rms_full
