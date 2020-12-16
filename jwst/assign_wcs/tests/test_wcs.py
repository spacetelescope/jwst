import pytest
from numpy import zeros
from numpy.testing import assert_allclose
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from asdf.tests import helpers

from jwst.assign_wcs import AssignWcsStep, pointing
from jwst.transforms import models
from jwst.datamodels import ImageModel, CubeModel, open
from gwcs.wcstools import grid_from_bounding_box


def create_hdul(wcskeys={
                      'wcsaxes': 2,
                      'ra_ref': 22.02351763251896,
                      'dec_ref': 11.99875540218638,
                      'v2_ref': 86.039011,
                      'v3_ref': -493.385704,
                      'roll_ref': 0.005076934167039675},
                data_shape=(2048, 2048)):

    """Create nircam hdulist specifically for the SIP test."""
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['DATAMODL'] = 'ImageModel'
    phdu.header['TELESCOP'] = "JWST"
    phdu.header['FILENAME'] = "test+F444W"
    phdu.header['INSTRUME'] = 'NIRCAM'
    phdu.header['CHANNEL'] = 'LONG'
    phdu.header['DETECTOR'] = 'NRCALONG'
    phdu.header['FILTER'] = 'F444W'
    phdu.header['PUPIL'] = 'CLEAR'
    phdu.header['MODULE'] = 'A'
    phdu.header['TIME-OBS'] = '16:58:27.258'
    phdu.header['DATE-OBS'] = '2021-10-25'
    phdu.header['EXP_TYPE'] = 'NRC_IMAGE'
    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header['SUBARRAY'] = 'FULL'

    scihdu.header.update(wcskeys)

    scihdu.data = zeros(data_shape)

    hdul.append(phdu)
    hdul.append(scihdu)

    return hdul


@pytest.fixture
def create_model_2d():
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


@pytest.fixture
def create_model_3d():
    im = CubeModel()
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
    assert_allclose(roll_angle, r0, atol=1e-3)


def test_v23_to_sky():
    """
    Test taken from INS report.
    """
    ra_ref = 165  # in deg
    dec_ref = 54  # in deg
    v2_ref = -503.654472 / 3600  # in deg
    v3_ref = -318.742464 / 3600  # in deg
    r0 = 37  # in deg

    v2 = 210  # in deg
    v3 = -75  # in deg
    expected_ra_dec = (107.12810484789563, -35.97940247128502)  # in deg
    angles = [-v2_ref, v3_ref, -r0, -dec_ref, ra_ref]
    axes = "zyxyz"
    v2s = models.V23ToSky(angles, axes_order=axes)
    radec = v2s(v2, v3)
    assert_allclose(radec, expected_ra_dec, atol=1e-10)


def test_frame_from_model_3d(tmpdir, create_model_3d):
    """ Tests creating a frame from a data model. """
    # Test CompositeFrame initialization (celestial and spectral)
    im = create_model_3d
    frame = pointing.frame_from_model(im)
    radec, lam = frame.coordinates(1, 2, 3)

    assert_allclose(radec.spherical.lon.value, 1)
    assert_allclose(radec.spherical.lat.value, 2)
    u.allclose(lam, 3 * u.um)

    # Test CompositeFrame initialization with custom frames
    im.meta.wcsinfo.ctype1 = 'ALPHA1A'
    im.meta.wcsinfo.ctype2 = 'BETA1A'
    frame = pointing.frame_from_model(im)

    assert frame.frames[1].name == 'ALPHA1A_BETA1A'
    assert frame.frames[1].axes_names == ('ALPHA1A', 'BETA1A')

    tree = {'frame': frame}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_frame_from_model_2d(tmpdir, create_model_2d):
    """ Tests creating a frame from a data model. """
    # Test 2D spatial custom frame
    im = create_model_2d
    frame = pointing.frame_from_model(im)

    assert frame.name == "sky"
    assert frame.axes_names == ("RA", "DEC")

    tree = {'frame': frame}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_create_fitswcs(tmpdir, create_model_3d):
    """GWCS from create_fitswcs function and astropy.wcs give same result"""
    im = create_model_3d
    w3d = pointing.create_fitswcs(im)
    gra, gdec, glam = w3d(1, 1, 1)

    path = str(tmpdir.join("fitswcs.fits"))
    im.save(path)
    with fits.open(path) as hdulist:
        hdu = hdulist["SCI"]
        w = wcs.WCS(hdu.header)
    wcel = w.sub(['celestial'])
    ra, dec = wcel.all_pix2world(1, 1, 0)

    # Check that astropy.wcs.WCS and gwcs.WCS give same result
    assert_allclose((ra, dec), (gra, gdec))

    # test serialization
    tree = {'wcs': w3d}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_sip_approx(tmpdir):
    # some of the wcs info
    true_wcs = {
                'ctype1': 'RA---TAN-SIP',
                'ctype2': 'DEC--TAN-SIP',
                'crpix1': 1025.0,
                'crpix2': 1025.0,
                'crval1': 22.023508708126627,
                'crval2': 11.998764139672561,
                'cd1_1': -1.7436649443640038e-05,
                'cd1_2': -2.0849272804435688e-08,
                'cd2_1': -4.616114099568815e-08,
                'cd2_2': 1.7519905482510155e-05,
                'a_order': 3,
                'b_order': 3,
                'a_0_2': -1.5520006349224245e-06,
                'a_0_3': 4.3714933936267264e-13,
                'a_1_1': -1.1604291350572094e-05,
                'a_1_2': 1.7125037791723067e-09,
                'a_2_0': 1.9225112320436695e-06,
                'a_2_1': -7.926931183344235e-11,
                'a_3_0': 1.5974716218080076e-09,
                'b_0_2': -6.7554035599749904e-06,
                'b_0_3': 1.6911444547233376e-09,
                'b_1_1': 3.6345928479129676e-06,
                'b_1_2': -8.64927324063469e-11,
                'b_2_0': 4.915411904251364e-06,
                'b_2_1': 1.5701011999755348e-09,
                'b_3_0': -2.9577497324389284e-12,
                'sipmxerr': 3.950294659161568e-11,
                'sipiverr': 0.11047793758675752,
                'ap_0_1': 3.988114493629691e-07,
                'ap_0_2': 1.5486628623564042e-06,
                'ap_0_3': 3.784084622814559e-11,
                'ap_1_0': -7.65706188115184e-06,
                'ap_1_1': 1.149989280312782e-05,
                'ap_1_2': -1.4994544929498265e-09,
                'ap_2_0': -1.904283418222981e-06,
                'ap_2_1': -4.507459210843236e-11,
                'ap_3_0': -1.630572028632995e-09,
                'bp_0_1': -7.46518175488653e-06,
                'bp_0_2': 6.690996263911794e-06,
                'bp_0_3': -1.5898038512674783e-09,
                'bp_1_0': 4.524878821020556e-07,
                'bp_1_1': -3.60126568682908e-06,
                'bp_1_2': -4.4691036548847474e-11,
                'bp_2_0': -4.9044713849991446e-06,
                'bp_2_1': -1.7110976247273783e-09,
                'bp_3_0': 3.895381500452847e-11,
                'ap_order': 3,
                'bp_order': 3}

    hdu1 = create_hdul()
    im = ImageModel(hdu1)

    pipe = AssignWcsStep()
    result = pipe.call(im)

    # check that result.meta.wcsinfo has correct
    # values after SIP approx.
    wcs_info = result.meta.wcsinfo.instance

    # make sure all expected keys are there
    assert set(true_wcs.keys()).issubset(set((wcs_info.keys())))

    # and make sure they match
    for key in true_wcs:
        if 'ctype' in key:
            assert true_wcs[key] == wcs_info[key]
        else:
            assert_allclose(true_wcs[key], wcs_info[key])

    # evaulate fits wcs and gwcs, make sure they agree
    grid = grid_from_bounding_box(result.meta.wcs.bounding_box)
    gwcs_ra, gwcs_dec = result.meta.wcs(*grid)
    fits_wcs = wcs.WCS(wcs_info)
    fitswcs_res = fits_wcs.pixel_to_world(*grid)

    assert_allclose(fitswcs_res.ra.deg, gwcs_ra)
    assert_allclose(fitswcs_res.dec.deg, gwcs_dec, atol=1.5e-6)

    # now write the file out, read it back in, and check that the fit values are preserved
    path = str(tmpdir.join("tmp_sip_wcs.fits"))
    result.write(path)


    with open(path) as result_read:
        result.meta.wcsinfo == result_read.meta.wcsinfo
