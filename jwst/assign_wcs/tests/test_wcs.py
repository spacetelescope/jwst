import pytest
from numpy.testing import utils
from astropy import units as u
from astropy import wcs
from astropy.tests.helper import  assert_quantity_allclose
from asdf.tests import helpers

from .. import pointing
from ...transforms import models
from ...datamodels import ImageModel, fits_support



def _create_model_2d():
    im = ImageModel()
    im.meta.wcsinfo.crpix1 = 2.5
    im.meta.wcsinfo.crpix2 = 3
    im.meta.wcsinfo.crval1 = 5.6
    im.meta.wcsinfo.crval2 = -72.3
    im.meta.wcsinfo.wcsaxes = 2
    im.meta.wcsinfo.cunit1 = 'deg'
    im.meta.wcsinfo.cunit2 = 'deg'
    im.meta.wcsinfo.ctype1 = 'RA---TAN'
    im.meta.wcsinfo.ctype2 = 'DEC--TAN'
    im.meta.wcsinfo.pc1_1 = 1.
    im.meta.wcsinfo.pc1_2 = 0
    im.meta.wcsinfo.pc2_1 = 0.
    im.meta.wcsinfo.pc2_2 = 1.
    return im


def _create_model_3d():
    im = _create_model_2d()
    im.meta.wcsinfo.crpix3 = 1
    im.meta.wcsinfo.crval3 = 100
    im.meta.wcsinfo.wcsaxes = 3
    im.meta.wcsinfo.cunit3 = 'um'
    im.meta.wcsinfo.ctype3 = 'WAVE'
    im.meta.wcsinfo.pc3_1 = 0.
    im.meta.wcsinfo.pc3_2 = 0
    im.meta.wcsinfo.pc3_3 = 1
    im.meta.wcsinfo.pc1_3 = 0
    im.meta.wcsinfo.pc2_3 = 0.
    return im


def test_roll_angle():
    """
    A sanity test - when V2_REF = 0 and V3_REF = 0,
    for V2, V3 close to he origin
    ROLL_REF should be approximately PA_V3 .

    (Test taken from SIAF report.)
    """
    ra_ref = 165 # in deg
    dec_ref = 54 # in deg
    v2_ref = 0
    v3_ref = 0
    r0 = 37 # in deg

    v2 = .01 # in arcsec
    v3 = .01 # in arcsec
    roll_angle = pointing.compute_roll_ref(v2_ref, v3_ref, r0, ra_ref, dec_ref, v2, v3)
    utils.assert_allclose(roll_angle, r0, atol=1e-3)


def test_v23_to_sky():
    """
    Test taken from INS report.
    """
    ra_ref = 165 # in deg
    dec_ref = 54 # in deg
    v2_ref = -503.654472 / 3600 # in deg
    v3_ref = -318.742464 / 3600 # in deg
    r0 = 37 # in deg

    v2 = 210 # in deg
    v3 = -75 # in deg
    expected_ra_dec = (107.12810484789563, -35.97940247128502) # in deg
    angles = [-v2_ref, v3_ref, -r0, -dec_ref, ra_ref]
    axes = "zyxyz"
    v2s = models.V23ToSky(angles, axes_order=axes)
    radec = v2s(v2, v3)
    utils.assert_allclose(radec, expected_ra_dec, atol=1e-10)


def test_frame_from_model(tmpdir):
    """ Tests creating a frame from a data model. """
    # Test CompositeFrame initialization (celestial and spectral)
    im = _create_model_3d()
    frame = pointing.frame_from_model(im)
    radec, lam = frame.coordinates(1, 2, 3)
    utils.assert_allclose(radec.spherical.lon.value, 1)
    utils.assert_allclose(radec.spherical.lat.value, 2)
    assert_quantity_allclose(lam, 3 * u.um)

    # Test CompositeFrame initialization with custom frames
    im.meta.wcsinfo.ctype1 = 'ALPHA1A'
    im.meta.wcsinfo.ctype2 = 'BETA1A'
    frame = pointing.frame_from_model(im)
    assert frame.frames[1].name == 'ALPHA1A_BETA1A'
    assert frame.frames[1].axes_names == ('ALPHA1A', 'BETA1A')

    tree = {'frame': frame}
    helpers.assert_roundtrip_tree(tree, tmpdir)

    # Test 2D spatial custom frame
    im = _create_model_2d()
    frame = pointing.frame_from_model(im)
    assert frame.name == "sky"
    assert frame.axes_names == ("RA", "DEC")
    tree = {'frame': frame}
    helpers.assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.filterwarnings("ignore: The WCS transformation has more axes")
def test_create_fitwcs(tmpdir):
    """ Test GWCS vs FITS WCS results. """
    im = _create_model_3d()
    w3d = pointing.create_fitswcs(im)
    gra, gdec, glam = w3d(1, 1, 1)

    ff = fits_support.to_fits(im._instance, im._schema)
    hdu = fits_support.get_hdu(ff._hdulist, "SCI")
    w = wcs.WCS(hdu.header)
    wcel = w.sub(['celestial'])
    ra, dec = wcel.all_pix2world(1, 1, 1)
    utils.assert_allclose((ra, dec), (gra, gdec))

    # test serialization
    tree = {'wcs': w3d}
    helpers.assert_roundtrip_tree(tree, tmpdir)
