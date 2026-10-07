"""
Test NIRCAM grism WCS transformations.

Notes:
These test the stability of the WFSS and TSO transformations based on the
the reference file that is returned from CRDS. The absolute validity
of the results is verified by the team based on the specwcs reference
file.

"""

import numpy as np
import pytest
from stcal.alignment.util import sregion_to_footprint

from jwst.assign_wcs import nircam, util
from jwst.tests import nircam_rate_helpers as helpers
from jwst.tests.wcs_helpers import get_reference_files

# Allowed settings for nircam
tsgrism_filters = ["F277W", "F444W", "F322W2", "F356W"]

nircam_wfss_frames = ["grism_detector", "detector", "v2v3", "v2v3vacorr", "world"]

nircam_tsgrism_frames = ["grism_detector", "direct_image", "v2v3", "v2v3vacorr", "world"]

nircam_imaging_frames = ["detector", "v2v3", "v2v3vacorr", "world"]


@pytest.fixture(scope="module", params=tsgrism_filters)
def tsgrism_inputs(request):
    def _add_missing_key(missing_key=None, missing_offsets=False):
        model = helpers.nircam_tsgrism_rate_model(filter_name=request.param, with_wcs=False)
        if missing_key is not None:
            setattr(model.meta.wcsinfo, missing_key, None)
        if missing_offsets:
            model.meta.dither.x_offset = None
            model.meta.dither.y_offset = None
        return model, get_reference_files(model)

    return _add_missing_key


def test_nircam_wfss_available_frames():
    """Make sure that the expected GWCS reference frames are created."""
    for p in ["GRISMR", "GRISMC"]:
        wcsobj = helpers.nircam_wfss_rate_model(p).meta.wcs
        available_frames = wcsobj.available_frames
        assert all([a == b for a, b in zip(nircam_wfss_frames, available_frames)])


def test_nircam_tso_available_frames():
    """Make sure that the expected GWCS reference frames for TSO are created."""
    wcsobj = helpers.nircam_tsgrism_rate_model().meta.wcs
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(nircam_tsgrism_frames, available_frames)])


@pytest.mark.parametrize("key", ["siaf_xref_sci", "siaf_yref_sci"])
def test_extract_tso_object_fails_without_xref_yref(tsgrism_inputs, key):
    """
    TypeError for xref_sci because computing dither offset is attempted, get float - NoneType
    ValueError for yref_sci because computing dither offset is not attempted
    """
    with pytest.raises((TypeError, ValueError)):
        nircam.tsgrism(*tsgrism_inputs(missing_key=key))


def traverse_wfss_trace(pupil):
    """Make sure that the WFSS dispersion polynomials are reversible."""
    wcsobj = helpers.nircam_wfss_rate_model(pupil).meta.wcs
    detector_to_grism = wcsobj.get_transform("detector", "grism_detector")
    grism_to_detector = wcsobj.get_transform("grism_detector", "detector")

    # check the round trip, grism pixel 100,100, source at 110,110,order 1
    xgrism, ygrism, xsource, ysource, order_in = (100, 100, 110, 110, 1)
    x0, y0, lam, order = grism_to_detector(xgrism, ygrism, xsource, ysource, order_in)
    x, y, xdet, ydet, orderdet = detector_to_grism(x0, y0, lam, order)

    assert x0 == xsource
    assert y0 == ysource
    assert order == order_in
    assert xdet == xsource
    assert ydet == ysource
    assert orderdet == order_in


def test_traverse_wfss_grisms():
    """Check the dispersion polynomials for each grism."""
    for pupil in ["GRISMR", "GRISMC"]:
        traverse_wfss_trace(pupil)


def test_traverse_tso_grism():
    """Make sure that the TSO dispersion polynomials are reversible.
    All assert statements are in pixel space so 1/1000 px seems easily acceptable"""
    wcsobj = helpers.nircam_tsgrism_rate_model().meta.wcs
    detector_to_grism = wcsobj.get_transform("direct_image", "grism_detector")
    grism_to_detector = wcsobj.get_transform("grism_detector", "direct_image")

    # TSGRISM always has same source locations
    # takes x,y,order -> ra, dec, wave, order
    xin, yin, order = (1024, 34.448, 1)  # y must be on the trace to round trip
    x0, y0, lam, orderdet = grism_to_detector(xin, yin, order)
    x, y, orderdet = detector_to_grism(x0, y0, lam, order)

    # Check returned reference positions
    assert np.isclose(x0, helpers.TSO_WCS_KW["xref_sci"] - 1, atol=0.5)
    # y reference includes a ~22.5 pixel shift from TA to science position,
    # from the dither offset
    assert np.isclose(y0, helpers.TSO_WCS_KW["yref_sci"] - 1 + 22.5, atol=0.5)
    assert order == orderdet

    # Check round trip
    assert np.isclose(x, xin, atol=2e-2)
    assert np.isclose(y, yin, atol=2e-2)


def test_tsgrism_offset_warning(caplog, tsgrism_inputs):
    # Run the pipeline on data missing x/y offsets
    nircam.tsgrism(*tsgrism_inputs(missing_offsets=True))

    # A warning is logged
    assert "could not be applied" in caplog.text
    assert "may be inaccurate" in caplog.text


def test_imaging_frames():
    """Verify the available imaging mode reference frames."""
    wcsobj = helpers.nircam_image_rate_model().meta.wcs
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(nircam_imaging_frames, available_frames)])


def test_wfss_sip():
    wfss_model = helpers.nircam_wfss_rate_model()
    util.wfss_imaging_wcs(wfss_model, nircam.imaging, bbox=((1, 1024), (1, 1024)))
    for key in ["a_order", "b_order", "crpix1", "crpix2", "crval1", "crval2", "cd1_1"]:
        assert key in wfss_model.meta.wcsinfo.instance


def test_update_s_region_imaging():
    """Ensure the s_region keyword matches output of wcs.footprint()"""
    model = helpers.nircam_image_rate_model()
    util.update_s_region_imaging(model)

    s_region = model.meta.wcsinfo.s_region
    footprint = model.meta.wcs.footprint().flatten()
    footprint_sregion = sregion_to_footprint(s_region).flatten()
    assert np.allclose(footprint, footprint_sregion, rtol=1e-9)
