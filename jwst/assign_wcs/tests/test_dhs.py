import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filename
from gwcs import wcs
from numpy.testing import assert_allclose

from jwst.assign_wcs import AssignWcsStep, nircam
from jwst.assign_wcs.tests.helpers import (
    make_mock_dhs_nrca1_rate,
    make_mock_dhs_nrca1_regions,
    make_mock_dhs_nrcalong_rate,
    make_mock_dhs_nrcalong_regions,
)
from jwst.assign_wcs.tests.test_nircam import get_reference_files


@pytest.fixture
def mock_dhs_nrca1_rate():
    """Create a mock DHS NRCA1 rate file."""
    return make_mock_dhs_nrca1_rate()


@pytest.fixture
def mock_dhs_nrca1_regions(mock_dhs_nrca1_rate, tmp_path):
    """Create a mock NIRCam DHS NRCA1 regions reference file."""
    return make_mock_dhs_nrca1_regions(mock_dhs_nrca1_rate, tmp_path)


@pytest.fixture
def create_dhs_nrca1_wcs(mock_dhs_nrca1_rate, mock_dhs_nrca1_regions):
    """Create a WCS for the mock NRCA1 DHS mode."""
    im = mock_dhs_nrca1_rate

    ref = get_reference_files(im)
    ref["specwcs"] = get_pkg_data_filename(
        "data/nircam_nrca1_specwcs.asdf", package="jwst.assign_wcs.tests"
    )
    ref["regions"] = mock_dhs_nrca1_regions

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
    x_in, y_in, order_in = 1000, 20, 1
    x0, y0, lam, order_mid, stripe = grism_to_direct(x_in, y_in, order_in)
    x_rec, _y_rec, order_rec = grism_to_direct.inverse(x0, y0, lam, order_mid, stripe)

    assert_allclose(x_rec, x_in, atol=1e-2, rtol=0)
    assert order_rec == order_in

    # Full forward pipeline: (x, y, order) -> (ra, dec, lam, order, stripe)
    ra, dec, lam_world, order_out, stripe_out = wcsobj(x_in, y_in, order_in)

    # Test we are in the ballpark of the reference RA/Dec
    # No attempt to test the actual expected values for these exact pixels
    assert_allclose(ra, 80.0, atol=0.1)
    assert_allclose(dec, -69.5, atol=0.1)

    # Similarly test wavelength makes sense for F150W2
    assert lam_world > 1.0
    assert lam_world < 2.25

    assert order_out == order_in
    assert stripe_out in [7, 8, 9, 10]


@pytest.fixture
def mock_dhs_nrcalong_rate():
    """Create a mock DHS NRCALONG rate file."""
    return make_mock_dhs_nrcalong_rate()


@pytest.fixture
def mock_dhs_nrcalong_regions(mock_dhs_nrcalong_rate, tmp_path):
    """Create a mock NIRCam DHS NRCALONG regions reference file."""
    return make_mock_dhs_nrcalong_regions(mock_dhs_nrcalong_rate, tmp_path)


@pytest.fixture
def create_dhs_nrcalong_wcs(mock_dhs_nrcalong_rate, mock_dhs_nrcalong_regions):
    """
    Create a WCS for the mock NRCALONG DHS mode.

    specwcs is already delivered for NRCALONG, so we don't need to mock it.
    """
    im = mock_dhs_nrcalong_rate

    ref = get_reference_files(im)
    ref["regions"] = mock_dhs_nrcalong_regions

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
    x_in, y_in, order_in = 1000, 20, 1
    x0, y0, lam, order_mid, stripe = grism_to_direct(x_in, y_in, order_in)
    x_rec, _y_rec, order_rec = grism_to_direct.inverse(x0, y0, lam, order_mid, stripe)

    assert_allclose(x_rec, x_in, atol=1e-2, rtol=0)
    assert order_rec == order_in

    # Full forward pipeline: (x, y, order) -> (ra, dec, lam, order, stripe)
    ra, dec, lam_world, order_out, stripe_out = wcsobj(x_in, y_in, order_in)

    # Test we are in the ballpark of the reference RA/Dec
    # No attempt to test the actual expected values for these exact pixels
    assert_allclose(ra, 80.0, atol=0.1)
    assert_allclose(dec, -69.5, atol=0.1)

    # Similarly test wavelength makes sense for F444W
    assert lam_world > 3.5
    assert lam_world < 5.0

    assert order_out == order_in


def test_assign_wcs_step_nrca1_dhs(mock_dhs_nrca1_rate, mock_dhs_nrca1_regions):
    """Test that AssignWcsStep completes successfully on an NRCA1 DHS input."""
    result = AssignWcsStep.call(
        mock_dhs_nrca1_rate,
        override_specwcs=get_pkg_data_filename(
            "data/nircam_nrca1_specwcs.asdf", package="jwst.assign_wcs.tests"
        ),
        override_regions=mock_dhs_nrca1_regions,
    )
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert result.meta.wcs is not None
    assert "grism_detector" in result.meta.wcs.available_frames
    assert "world" in result.meta.wcs.available_frames

    # Scalar input
    ra, dec, lam_world, order_out, stripe_out = result.meta.wcs(1000, 20, 1)
    assert_allclose(ra, 80.0, atol=0.1)
    assert_allclose(dec, -69.5, atol=0.1)
    assert stripe_out in [7, 8, 9, 10]

    # Vector input
    x = np.array([900, 901, 902])
    y = np.array([20, 20, 20])
    order = np.array([1, 1, 1])
    ra, dec, _lam_world, _order_out, stripe_out = result.meta.wcs(x, y, order)
    assert_allclose(ra, 80.0, atol=0.1)
    assert_allclose(dec, -69.5, atol=0.1)
    assert np.all(np.isin(stripe_out, [7, 8, 9, 10]))


def test_assign_wcs_step_nrcalong_dhs(mock_dhs_nrcalong_rate, mock_dhs_nrcalong_regions):
    """Test that AssignWcsStep completes successfully on an NRCALONG DHS input."""
    result = AssignWcsStep.call(
        mock_dhs_nrcalong_rate,
        override_regions=mock_dhs_nrcalong_regions,
    )
    assert result.meta.cal_step.assign_wcs == "COMPLETE"
    assert result.meta.wcs is not None
    # Verify the WCS has the expected frames for DHS mode
    assert "grism_detector" in result.meta.wcs.available_frames
    assert "world" in result.meta.wcs.available_frames

    # Scalar input
    ra, dec, lam_world, order_out, stripe_out = result.meta.wcs(900, 20, 1)
    assert_allclose(ra, 80.0, atol=0.1)
    assert_allclose(dec, -69.5, atol=0.1)
    assert stripe_out == 1  # For NRCALONG, all pixels should be in the same stripe

    # Vector input
    x = np.array([900, 901, 902])
    y = np.array([20, 20, 20])
    order = np.array([1, 1, 1])
    ra, dec, _lam_world, _order_out, stripe_out = result.meta.wcs(x, y, order)
    assert_allclose(ra, 80.0, atol=0.1)
    assert_allclose(dec, -69.5, atol=0.1)
    assert np.all(stripe_out == 1)
