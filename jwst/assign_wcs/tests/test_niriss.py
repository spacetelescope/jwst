"""
Test NIRISS grism WCS transformations.

Notes:
These test the stability of the WFSS transformations based on the
the reference file that is returned from CRDS. The absolute validity
of the results is verified by the team based on the specwcs reference
file.

"""

from numpy.testing import assert_allclose

from jwst.assign_wcs import niriss, util
from jwst.tests import niriss_rate_helpers as helpers

# Allowed settings for niriss
niriss_wfss_frames = ["grism_detector", "detector", "v2v3", "v2v3vacorr", "world"]
niriss_imaging_frames = ["detector", "v2v3", "v2v3vacorr", "world"]
niriss_grisms = ["GR150R", "GR150C"]


def test_niriss_wfss_available_frames():
    for f in ["GR150R", "GR150C"]:
        wcsobj = helpers.niriss_wfss_rate_model(f).meta.wcs
        available_frames = wcsobj.available_frames
        assert all([a == b for a, b in zip(niriss_wfss_frames, available_frames)])


def traverse_wfss_trace(filter_name):
    wcsobj = helpers.niriss_wfss_rate_model(filter_name=filter_name).meta.wcs
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
    """Make sure the trace polynomials roundtrip for both grisms."""
    for f in niriss_grisms:
        traverse_wfss_trace(f)


# 1.585 is 10x the angular repeatability of the niriss pupil wheel
def test_filter_rotation(theta=[-0.1, 0, 0.5, 1.585]):
    """Make sure that the filter rotation is reversible."""
    for f in niriss_grisms:
        wcsobj = helpers.niriss_wfss_rate_model(f).meta.wcs
        g2d = wcsobj.get_transform("grism_detector", "detector")
        d2g = wcsobj.get_transform("detector", "grism_detector")
        for angle in theta:
            d2g.theta = angle
            g2d.theta = -angle
            xsource, ysource, wave, order = (110, 110, 2.3, 1)
            xgrism, ygrism, xs, ys, orderout = d2g(xsource, ysource, wave, order)
            xsdet, ysdet, wavedet, orderdet = g2d(xgrism, ygrism, xs, ys, orderout)
            assert_allclose(xsdet, xsource)
            assert_allclose(ysdet, ysource)
            assert_allclose(wavedet, wave, atol=2e-3)
            assert orderdet == order


def test_imaging_frames():
    """Verify the available imaging mode reference frames."""
    wcsobj = helpers.niriss_image_rate_model().meta.wcs
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(niriss_imaging_frames, available_frames)])


def test_imaging_distortion():
    """Verify that the distortion correction roundtrips."""
    wcsobj = helpers.niriss_image_rate_model().meta.wcs
    sky_to_detector = wcsobj.get_transform("world", "detector")
    detector_to_sky = wcsobj.get_transform("detector", "world")

    # we'll use the crpix as the simplest reference point
    ra = helpers.DEFAULT_WCS_KW["crval1"]
    dec = helpers.DEFAULT_WCS_KW["crval2"]

    x, y = sky_to_detector(ra, dec)
    raout, decout = detector_to_sky(x, y)

    assert_allclose(raout, ra)
    assert_allclose(decout, dec)


def test_wfss_sip():
    wfss_model = helpers.niriss_wfss_rate_model()
    util.wfss_imaging_wcs(
        wfss_model, niriss.imaging, max_pix_error=0.05, bbox=((1, 1024), (1, 1024))
    )
    for key in ["a_order", "b_order", "crpix1", "crpix2", "crval1", "crval2", "cd1_1"]:
        assert key in wfss_model.meta.wcsinfo.instance
