"""Unit tests for AMI bp_fix module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
import stdatamodels.jwst.datamodels as dm

from jwst.ami import bp_fix, utils


@pytest.mark.parametrize("filt", bp_fix.filts)
def test_create_wavelengths(filt):

    wl_ctr, wl_left, wl_right = bp_fix.create_wavelengths(filt)

    expected_ctr = bp_fix.filtwl_d[filt]
    halfpower_left = bp_fix.filthp_d[filt][0]
    halfpower_right = bp_fix.filthp_d[filt][1]
    assert np.isclose(wl_ctr, expected_ctr)
    assert wl_left < halfpower_left
    assert wl_right > halfpower_right


@pytest.fixture
def circular_pupil():
    shape = (1024, 1024)
    r = 0.2
    x = np.linspace(-1, 1, shape[0])
    y = np.linspace(-1, 1, shape[1])
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    pupil = np.zeros(shape)
    pupil[rr < r] = 1
    return pupil


@pytest.fixture
def nrm_model(circular_pupil):

    return dm.NRMModel(nrm=circular_pupil)


def test_calc_psf(example_model, circular_pupil):
    """
    Tests for calc_psf and calc_pupil_support.
    
    These have very similar inputs and outputs so are tested together.
    """

    filt = example_model.meta.instrument.filter
    wl = bp_fix.filtwl_d[filt]
    pxsc_deg = utils.degrees_per_pixel(example_model)[0]
    pxsc_rad = pxsc_deg * np.pi / (180)
    fov_npix = example_model.data.shape[1]

    im_single_psf = bp_fix.calcpsf(wl, fov_npix, pxsc_rad, circular_pupil)
    im_pupil = bp_fix.calc_pupil_support(filt, fov_npix, pxsc_rad, circular_pupil)

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


@pytest.fixture
def data_with_bad_pixels(example_model):
    rng = np.random.default_rng(0)
    shp = example_model.data.shape
    data = rng.normal(size=shp)
    data[0, 20, 20] = 100
    return data


def test_fourier_corr(data_with_bad_pixels):
    """Test bad_pixels and fourier_corr functions."""

    data = data_with_bad_pixels[0]
    badpx = bp_fix.bad_pixels(data, 8, 10)
    assert badpx.shape == data.shape
    assert badpx[20, 20] == 1
    assert np.sum(badpx) == 1
    assert badpx.dtype == np.int64

    # fmas needs to be half-size
    fmas = bp_fix.transform_image(data)[:, : data.shape[1] // 2 + 1]
    data_fixed = bp_fix.fourier_corr(data, badpx, fmas)

    assert data_fixed.shape == data.shape
    assert not np.isnan(data_fixed).any()
    assert not np.isinf(data_fixed).any()
    assert data_fixed.dtype == np.float64

    # check that the bad pixel was fixed. should be within 1 sigma of the noise
    # since we are using a relatively large median region
    assert np.abs(data_fixed[20, 20]) < 1


def test_fix_bad_pixels(example_model, data_with_bad_pixels, nrm_model):

    data = data_with_bad_pixels
    pxdq0 = np.zeros_like(data, dtype=np.bool_)
    filt = example_model.meta.instrument.filter  
    pxsc_deg = utils.degrees_per_pixel(example_model)[0]
    pxsc_mas = pxsc_deg * 3600000

    data_out, pxdq_out = bp_fix.fix_bad_pixels(data, pxdq0, filt, pxsc_mas, nrm_model)

    assert data_out.shape == data.shape
    assert pxdq_out.shape == data.shape
    assert pxdq_out[0, 20, 20] == 1
    assert np.sum(pxdq_out) == 1
    assert data_out.dtype == np.float64
    assert pxdq_out.dtype == np.int64
    
    # replace bad pixel in both arrays with a placeholder so can assert the rest did not change
    data[0, 20, 20] = 0
    data_out[0, 20, 20] = 0
    assert_allclose(data_out, data, rtol=1e-5)
