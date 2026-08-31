import numpy as np
import pytest
from gwcs import wcs
from numpy.testing import assert_allclose

from jwst.assign_wcs import AssignWcsStep, nircam
from jwst.assign_wcs.tests.helpers import (
    make_mock_dhs_nrca1_rate,
    make_mock_dhs_nrcalong_rate,
)
from jwst.assign_wcs.tests.test_nircam import get_reference_files


@pytest.fixture
def mock_dhs_nrca1_rate():
    """Create a mock DHS NRCA1 rate file."""
    return make_mock_dhs_nrca1_rate()


@pytest.fixture
def create_dhs_nrca1_wcs(mock_dhs_nrca1_rate):
    """Create a WCS for the mock NRCA1 DHS mode."""
    im = mock_dhs_nrca1_rate

    ref = get_reference_files(im)
    pipeline = nircam.dhs(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def test_dhs_nrca1_roundtrip(create_dhs_nrca1_wcs):
    """
    Test that the DHS mode WCS round-trips.

    Also verifies that:
    - Sky coordinates from the full forward pipeline are physically near the
      reference pointing and in a plausible spectral range.
    - Spectral dispersion varies along the x axis.
    - The stripe ID is one of the expected DHS stripe IDs.
    """
    wcsobj = create_dhs_nrca1_wcs
    grism_to_direct = wcsobj.get_transform("grism_detector", "direct_image")

    # Roundtrip: grism_detector -> direct_image -> grism_detector
    # Only the x (dispersion) coordinate is expected to round-trip, because the
    # forward transform collapses all y positions to the source reference point.
    x_in = [1000] * 4  # same x for all stripes
    y_in = 32.5 + np.arange(0, 260, 65)  # center y for each stripe: (32.5, 97.5, 162.5, 227.5)
    order_in = [1] * 4  # same order for all stripes
    x0, y0, lam, order_mid, stripe = grism_to_direct(x_in, y_in, order_in)
    x_rec, y_rec, order_rec = grism_to_direct.inverse(x0, y0, lam, order_mid, stripe)

    assert_allclose(x_rec, x_in, atol=1e-2, rtol=0)
    assert_allclose(order_rec, order_in)

    # The y coordinate returned should be at the expected trace position for each stripe:
    # these are at the center of the stripe, but may be a few pixels up or down.
    # Current calibration is preliminary, so allow a couple pixels tolerance on the
    # expected value.
    expected = [34, 105, 171, 223]
    assert_allclose(y_rec, expected, atol=2)
    assert_allclose(y_rec, y_in, atol=10)

    # Full forward pipeline: (x, y, order) -> (ra, dec, lam, order, stripe)
    ra, dec, lam_world, order_out, stripe_out = wcsobj(x_in, y_in, order_in)

    # Test we are in the ballpark of the reference RA/Dec.
    # Values for all stripes are the same.
    assert_allclose(ra, 265.76, atol=0.02)
    assert_allclose(dec, 66.93, atol=0.02)

    # Similarly test wavelength makes sense for F150W2
    assert np.all(lam_world > 1.0)
    assert np.all(lam_world < 2.25)

    # Check expected order and stripe values
    assert_allclose(order_rec, order_in)
    assert_allclose(stripe_out, [10, 9, 8, 7])


@pytest.fixture
def mock_dhs_nrcalong_rate():
    """Create a mock DHS NRCALONG rate file."""
    return make_mock_dhs_nrcalong_rate()


@pytest.fixture
def create_dhs_nrcalong_wcs(mock_dhs_nrcalong_rate):
    """
    Create a WCS for the mock NRCALONG DHS mode.
    """
    im = mock_dhs_nrcalong_rate
    ref = get_reference_files(im)

    pipeline = nircam.dhs(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def test_dhs_nrcalong_roundtrip(create_dhs_nrcalong_wcs):
    """
    Test that the DHS mode WCS round-trips for NRCALONG.

    Also verifies that:
    - Sky coordinates from the full forward pipeline are physically near the
      reference pointing and in a plausible spectral range.
    - Spectral dispersion varies along the x axis.
    """
    wcsobj = create_dhs_nrcalong_wcs
    grism_to_direct = wcsobj.get_transform("grism_detector", "direct_image")

    # Roundtrip: grism_detector -> direct_image -> grism_detector
    # Only the x (dispersion) coordinate is expected to round-trip, because the
    # forward transform collapses all y positions to the source reference point.
    x_in = [1000] * 4  # same x for all stripes
    y_in = 32.5 + np.arange(0, 260, 65)  # center y for each stripe: (32.5, 97.5, 162.5, 227.5)
    order_in = [1] * 4  # same order for all stripes
    x0, y0, lam, order_mid, stripe = grism_to_direct(x_in, y_in, order_in)
    x_rec, y_rec, order_rec = grism_to_direct.inverse(x0, y0, lam, order_mid, stripe)

    assert_allclose(x_rec, x_in, atol=1e-2, rtol=0)
    assert_allclose(order_rec, order_in)

    # The y coordinate returned should be at the expected trace position for each stripe:
    # these are at the center of the stripe, so they match the input.
    assert_allclose(y_rec, y_in, atol=1e-2, rtol=0)

    # Full forward pipeline: (x, y, order) -> (ra, dec, lam, order, stripe)
    ra, dec, lam_world, order_out, stripe_out = wcsobj(x_in, y_in, order_in)

    # Test we are in the ballpark of the reference RA/Dec
    # Values for all stripes are the same.
    assert_allclose(ra, 265.76, atol=0.02)
    assert_allclose(dec, 66.93, atol=0.02)

    # Similarly test wavelength makes sense for F332W2
    assert np.all(lam_world > 2.0)
    assert np.all(lam_world < 4.5)

    # Check expected order and stripe values
    assert_allclose(order_rec, order_in)
    assert_allclose(stripe_out, [1, 2, 3, 4])


def test_assign_wcs_step_nrca1_dhs(mock_dhs_nrca1_rate):
    """Test that AssignWcsStep completes successfully on an NRCA1 DHS input."""
    result = AssignWcsStep.call(mock_dhs_nrca1_rate)
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert result.meta.wcs is not None
    assert "grism_detector" in result.meta.wcs.available_frames
    assert "world" in result.meta.wcs.available_frames

    # Scalar input
    ra, dec, lam_world, order_out, stripe_out = result.meta.wcs(1000, 20, 1)
    assert_allclose(ra, 265.7, atol=0.1)
    assert_allclose(dec, 66.9, atol=0.1)
    assert stripe_out in [7, 8, 9, 10]

    # Vector input
    x = np.array([900, 901, 902])
    y = np.array([20, 20, 20])
    order = np.array([1, 1, 1])
    ra, dec, _lam_world, _order_out, stripe_out = result.meta.wcs(x, y, order)
    assert_allclose(ra, 265.7, atol=0.1)
    assert_allclose(dec, 66.9, atol=0.1)
    assert np.all(np.isin(stripe_out, [7, 8, 9, 10]))


def test_assign_wcs_step_nrcalong_dhs(mock_dhs_nrcalong_rate):
    """Test that AssignWcsStep completes successfully on an NRCALONG DHS input."""
    result = AssignWcsStep.call(mock_dhs_nrcalong_rate)
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert result.meta.wcs is not None
    # Verify the WCS has the expected frames for DHS mode
    assert "grism_detector" in result.meta.wcs.available_frames
    assert "world" in result.meta.wcs.available_frames

    # Scalar input
    ra, dec, lam_world, order_out, stripe_out = result.meta.wcs(900, 20, 1)
    assert_allclose(ra, 265.7, atol=0.1)
    assert_allclose(dec, 66.9, atol=0.1)
    assert stripe_out == 1  # For NRCALONG, all pixels should be in the same stripe

    # Vector input
    x = np.array([900, 901, 902])
    y = np.array([20, 20, 20])
    order = np.array([1, 1, 1])
    ra, dec, _lam_world, _order_out, stripe_out = result.meta.wcs(x, y, order)
    assert_allclose(ra, 265.7, atol=0.1)
    assert_allclose(dec, 66.9, atol=0.1)
    assert np.all(stripe_out == 1)
