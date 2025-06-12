"""Unit tests for AMI bp_fix module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from jwst.ami import bp_fix
from .conftest import PXSC_RAD, PXSC_MAS


@pytest.mark.parametrize("filt", bp_fix.filts)
def test_create_wavelengths(filt):
    wl_ctr, wl_left, wl_right = bp_fix.create_wavelengths(filt)

    expected_ctr = bp_fix.filtwl_d[filt]
    halfpower_left = bp_fix.filthp_d[filt][0]
    halfpower_right = bp_fix.filthp_d[filt][1]
    assert np.isclose(wl_ctr, expected_ctr)
    assert wl_left < halfpower_left
    assert wl_right > halfpower_right


def test_calc_psf(example_model, circular_pupil):
    """
    Tests for calc_psf and calc_pupil_support.

    These have very similar inputs and outputs so are tested together.
    """

    filt = example_model.meta.instrument.filter
    wl = bp_fix.filtwl_d[filt]
    fov_npix = example_model.data.shape[1]

    im_single_psf = bp_fix.calcpsf(wl, fov_npix, PXSC_RAD, circular_pupil)
    im_pupil = bp_fix.calc_pupil_support(filt, fov_npix, PXSC_RAD, circular_pupil)

    for im in [im_single_psf, im_pupil]:
        # basic checks: no NaNs or infs, right shape and data type
        assert not np.any(np.isnan(im))
        assert not np.any(np.isinf(im))
        assert im.shape == (fov_npix, fov_npix)
        assert im.dtype == np.float64

        # output of calc_psf should be an Airy ring.
        # output of calc_pupil_support should be a fuzzy Airy ring.
        # These shapes are a bit annoying to check fully without recalculating them.
        # But can at least check that the peak is in the center and all positive
        assert np.min(im) >= 0
        assert np.max(im) == im[fov_npix // 2, fov_npix // 2]


def test_fourier_corr(example_model):
    """
    Test bad_pixels and fourier_corr functions.

    Example model has a bad pixel at (0, 20, 20).
    """

    data = example_model.data[0]
    badpx = bp_fix.bad_pixels(data, 3, 20)
    assert badpx.shape == data.shape
    assert badpx[20, 20] == 1
    assert np.sum(badpx) == 1
    assert badpx.dtype == np.bool_

    # fmas needs to be half-size
    fmas = bp_fix.transform_image(data)[:, : data.shape[1] // 2 + 1]
    data_fixed = bp_fix.fourier_corr(data, badpx, fmas)

    assert data_fixed.shape == data.shape
    assert not np.isnan(data_fixed).any()
    assert not np.isinf(data_fixed).any()
    assert data_fixed.dtype == np.float32

    # check that the bad pixel was fixed. should be within 1 sigma of the noise
    # since we are using a relatively large median region
    assert np.abs(data_fixed[20, 20]) < 1


def test_fix_bad_pixels(example_model, nrm_model_circular):
    """
    Test fix_bad_pixels function.

    Example model has a bad pixel at (0, 20, 20).
    """
    data = example_model.data.copy()
    nrm_model = nrm_model_circular
    pxdq0 = np.zeros_like(data, dtype=np.bool_)
    filt = example_model.meta.instrument.filter

    data_out, pxdq_out = bp_fix.fix_bad_pixels(data, pxdq0, filt, PXSC_MAS, nrm_model)

    assert data_out.shape == data.shape
    assert pxdq_out.shape == data.shape
    assert pxdq_out[0, 20, 20] == 1
    assert np.sum(pxdq_out) == 1
    assert data_out.dtype == np.float32
    assert pxdq_out.dtype == np.int64

    # replace bad pixel in both arrays with a placeholder so can assert the rest did not change
    data[0, 20, 20] = 0
    data_out[0, 20, 20] = 0
    assert_allclose(data_out, data, rtol=1e-5)


def test_fix_bad_pixels_even_size(example_model, nrm_model_circular):
    """Test coverage for bug where even-sized arrays were failing."""

    data = example_model.data.copy()
    nrm_model = nrm_model_circular
    pxdq0 = np.zeros_like(data, dtype=np.bool_)
    filt = example_model.meta.instrument.filter

    # Make the data even-sized
    data = data[:, : data.shape[1] - 1, : data.shape[2] - 1]
    pxdq0 = pxdq0[:, : pxdq0.shape[1] - 1, : pxdq0.shape[2] - 1]

    data_out, pxdq_out = bp_fix.fix_bad_pixels(data, pxdq0, filt, PXSC_MAS, nrm_model)

    assert data_out.shape == data.shape
    assert pxdq_out.shape == data.shape
    assert data_out[0, 20, 20] < 10
    assert pxdq_out[0, 20, 20] == 1
    assert np.sum(pxdq_out) == 1
    assert data_out.dtype == np.float32
    assert pxdq_out.dtype == np.int64

    # replace bad pixel in both arrays with a placeholder so can assert the rest did not change
    data[0, 20, 20] = 0
    data_out[0, 20, 20] = 0
    assert_allclose(data_out, data, rtol=1e-5)  # Use assert_allclose to compare arrays
