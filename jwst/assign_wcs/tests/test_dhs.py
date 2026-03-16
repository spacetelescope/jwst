import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.utils.data import get_pkg_data_filename
from gwcs import wcs
from numpy.testing import assert_allclose

from jwst.assign_wcs import nircam
from jwst.assign_wcs.tests.test_nircam import get_reference_files


@pytest.fixture
def mock_dhs_nrca1_rate():
    """Create a mock DHS NRCA1 rate file."""
    model = dm.CubeModel((5, 164, 2048))

    # Exposure
    model.meta.exposure.type = "NRC_TSGRISM"
    model.meta.exposure.ngroups = 5

    # Instrument
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "SHORT"
    model.meta.instrument.detector = "NRCA1"
    model.meta.instrument.filter = "F150W2"
    model.meta.instrument.pupil = "GDHS0"
    model.meta.instrument.module = "A"

    # Observation
    model.meta.observation.date = "2026-05-03"
    model.meta.observation.time = "00:00:00.000"

    # Subarray
    model.meta.subarray.name = "SUB164STRIPE4_DHS"
    model.meta.subarray.fastaxis = -1
    model.meta.subarray.slowaxis = 2
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 40
    model.meta.subarray.multistripe_skips1 = 1526
    model.meta.subarray.multistripe_skips2 = 85
    model.meta.subarray.num_superstripe = 0
    model.meta.subarray.repeat_stripe = 1
    model.meta.subarray.xstart = 1
    model.meta.subarray.xsize = 2048
    model.meta.subarray.ystart = 1
    model.meta.subarray.ysize = 164

    # WCS pointing
    model.meta.wcsinfo.ra_ref = 80.0
    model.meta.wcsinfo.dec_ref = -69.5
    model.meta.wcsinfo.v2_ref = 120.576793
    model.meta.wcsinfo.v3_ref = -527.501431
    model.meta.wcsinfo.roll_ref = 305.03225951982046
    model.meta.wcsinfo.velosys = -9378.83

    # Velocity aberration
    model.meta.velocity_aberration.scale_factor = 1.0

    return model


@pytest.fixture
def create_dhs_nrca1_wcs(mock_dhs_nrca1_rate):

    im = mock_dhs_nrca1_rate

    ref = get_reference_files(im)
    ref["specwcs"] = get_pkg_data_filename(
        "data/nircam_nrca1_specwcs.asdf", package="jwst.assign_wcs.tests"
    )
    ref["regions"] = get_pkg_data_filename(
        "data/tpauly_nrca1_regions.asdf", package="jwst.assign_wcs.tests"
    )

    pipeline = nircam.dhs(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def test_dhs_roundtrip(create_dhs_nrca1_wcs):
    """Test the DHS WCS dispersion roundtrip through intermediate frames.

    Transforms grism_detector → direct_image → grism_detector and checks that
    the x (dispersion-axis) pixel coordinate is recovered. This pattern mirrors
    ``test_traverse_tso_grism`` for the standard TSO mode.

    The y (cross-dispersion) pixel is not part of the roundtrip because the
    forward transform collapses all y positions to the source reference point
    (y0 = 0 for NRC_SHORT DHS), so y information is intentionally lost in the
    direct-image frame.

    Also verifies that:
    - Sky coordinates from the full forward pipeline are physically near the
      reference pointing and in a plausible NIRCam spectral range.
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
    x_rec, y_rec, order_rec = grism_to_direct.inverse(x0, y0, lam, order_mid, stripe)

    assert_allclose(x_rec, x_in, atol=1e-3)
    assert order_rec == order_in

    # Source should be pinned to the reference position (xc=0, yc=0 for NRC_SHORT DHS)
    assert_allclose(x0, 0.0, atol=1e-6)
    assert_allclose(y0, 0.0, atol=1e-6)

    # Full forward pipeline: (x, y, order) -> (ra, dec, lam, order, stripe)
    ra, dec, lam_world, order_out, stripe_out = wcsobj(x_in, y_in, order_in)
    assert np.isfinite(ra)
    assert np.isfinite(dec)
    assert np.isfinite(lam_world)
    assert_allclose(ra, 80.0, atol=1.0)
    assert_allclose(dec, -69.5, atol=1.0)
    assert 0.5 < lam_world < 5.0
    assert order_out == order_in
    assert stripe_out in [7, 8, 9, 10]

    # Spectral dispersion: different x positions yield different wavelengths
    _, _, lam_blue, _, _ = wcsobj(500, y_in, order_in)
    _, _, lam_red, _, _ = wcsobj(1500, y_in, order_in)
    assert lam_blue != lam_red
