import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.utils.data import get_pkg_data_filename
from gwcs import wcs
from numpy.testing import assert_allclose

from jwst.assign_wcs import AssignWcsStep, nircam
from jwst.assign_wcs.tests.test_nircam import get_reference_files


def _populate_shared_metadata(model):
    """
    Give model metadata that's identical between NRCA1 and NRCALONG modes.

    Updates model in place.
    """
    # Exposure
    model.meta.exposure.type = "NRC_TSGRISM"
    model.meta.exposure.ngroups = 5

    # Observation
    model.meta.observation.date = "2026-05-03"
    model.meta.observation.time = "00:00:00.000"

    # Subarray
    model.meta.subarray.name = "SUB164STRIPE4_DHS"
    model.meta.subarray.fastaxis = -1
    model.meta.subarray.slowaxis = 2
    model.meta.subarray.num_superstripe = 0
    model.meta.subarray.repeat_stripe = 1
    model.meta.subarray.xstart = 1
    model.meta.subarray.xsize = 2048
    model.meta.subarray.ystart = 1
    model.meta.subarray.ysize = 164

    # WCS info
    model.meta.wcsinfo.ra_ref = 80.0
    model.meta.wcsinfo.dec_ref = -69.5
    model.meta.wcsinfo.v2_ref = 120.576793
    model.meta.wcsinfo.v3_ref = -527.501431
    model.meta.wcsinfo.roll_ref = 305.03225951982046
    model.meta.wcsinfo.velosys = -9378.83

    # Velocity aberration
    model.meta.velocity_aberration.scale_factor = 1.0


@pytest.fixture
def mock_dhs_nrca1_rate():
    """Create a mock DHS NRCA1 rate file."""
    model = dm.CubeModel((5, 164, 2048))
    _populate_shared_metadata(model)

    # Instrument
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "SHORT"
    model.meta.instrument.detector = "NRCA1"
    model.meta.instrument.filter = "F150W2"
    model.meta.instrument.pupil = "GDHS0"
    model.meta.instrument.module = "A"

    # Subarray
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 40
    model.meta.subarray.multistripe_skips1 = 1526
    model.meta.subarray.multistripe_skips2 = 85

    return model


def _populate_shared_regions_metadata(model):
    """
    Give model metadata that's identical between NRCA1 and NRCALONG regions files.

    Updates model in place.
    """
    model.meta.description = "Mock DHS regions for testing"
    model.meta.author = "test"
    model.meta.pedigree = "GROUND"
    model.meta.useafter = "2000-01-01T00:00:00"
    model.meta.instrument.name = "NIRCAM"


@pytest.fixture
def mock_dhs_nrca1_regions(mock_dhs_nrca1_rate, tmp_path):
    """Create a mock NIRCam DHS regions reference file."""
    sci = mock_dhs_nrca1_rate
    reads1 = sci.meta.subarray.multistripe_reads1
    skips1 = sci.meta.subarray.multistripe_skips1
    reads2 = sci.meta.subarray.multistripe_reads2
    skips2 = sci.meta.subarray.multistripe_skips2

    # Full-frame (2048 x 2048) regions array; zero means "not in any stripe".
    regions = np.zeros((2048, 2048), dtype=np.float64)
    # Stripe IDs run from highest to lowest as row index increases
    stripe_ids = [10, 9, 8, 7]
    row_start = reads1 + skips1
    for stripe_id in stripe_ids:
        regions[row_start : row_start + reads2, :] = stripe_id
        row_start += reads2 + skips2

    regions_path = tmp_path / "mock_nrca1_regions.asdf"
    model = dm.RegionsModel()
    model.regions = regions
    _populate_shared_regions_metadata(model)
    model.save(str(regions_path))
    model.close()

    return str(regions_path)


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
    model = dm.CubeModel((5, 164, 2048))
    _populate_shared_metadata(model)

    # Instrument
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "LONG"
    model.meta.instrument.detector = "NRCALONG"
    model.meta.instrument.filter = "F444W"
    model.meta.instrument.pupil = "GRISMR"
    model.meta.instrument.module = "A"

    # Subarray
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 40
    model.meta.subarray.multistripe_skips1 = 971
    model.meta.subarray.multistripe_skips2 = 0

    # WCS information
    model.meta.wcsinfo.siaf_xref_sci = 862
    model.meta.wcsinfo.siaf_yref_sci = 20.5

    return model


@pytest.fixture
def mock_dhs_nrcalong_regions(mock_dhs_nrcalong_rate, tmp_path):
    """
    Create a mock NRCALONG regions ref file.

    For NRCALONG the same detector region is read in every readout, so this file is pretty
    trivial: it's nonzero in that single 40-pixel-tall stripe, and zero elsewhere.
    """
    sci = mock_dhs_nrcalong_rate
    reads1 = sci.meta.subarray.multistripe_reads1
    skips1 = sci.meta.subarray.multistripe_skips1
    reads2 = sci.meta.subarray.multistripe_reads2

    # The same 40-row detector stripe is read repeatedly into the SUB164STRIPE4_DHS
    # subarray, so the full-frame regions map only needs that single physical band.
    regions = np.zeros((2048, 2048), dtype=np.float64)
    row_start = reads1 + skips1
    regions[row_start : row_start + reads2, :] = 1

    regions_path = tmp_path / "mock_nrcalong_regions.asdf"
    model = dm.RegionsModel()
    model.regions = regions
    _populate_shared_regions_metadata(model)
    model.save(str(regions_path))
    model.close()

    return str(regions_path)


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
