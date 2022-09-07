import pytest

from numpy import zeros
from numpy.testing import assert_allclose
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.modeling.models import RotationSequence3D

from gwcs.wcstools import grid_from_bounding_box
from gwcs.geometry import SphericalToCartesian, CartesianToSpherical

from jwst.assign_wcs import AssignWcsStep, pointing
from jwst.datamodels import ImageModel, CubeModel, open


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
    for V2, V3 close to the origin
    ROLL_REF should be approximately PA_V3 .

    (Test taken from SIAF report.)
    """
    ra_ref = 165  # in deg
    dec_ref = 54  # in deg
    v2_ref = 0
    v3_ref = 0
    r0 = 37  # in deg

    v2 = .01  # in arcsec
    v3 = .01  # in arcsec
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
    angles = [v2_ref, -v3_ref, r0, dec_ref, -ra_ref]
    axes = "zyxyz"

    rot = RotationSequence3D(angles, axes_order=axes)
    v2s = SphericalToCartesian() | rot | CartesianToSpherical()
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


def test_frame_from_model_2d(tmpdir, create_model_2d):
    """ Tests creating a frame from a data model. """
    # Test 2D spatial custom frame
    im = create_model_2d
    frame = pointing.frame_from_model(im)

    assert frame.name == "sky"
    assert frame.axes_names == ("RA", "DEC")


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


def test_sip_approx(tmpdir):
    # some of the wcs info
    true_wcs = {
        'ctype1': 'RA---TAN-SIP',
        'ctype2': 'DEC--TAN-SIP',
        'crpix1': 1024.5,
        'crpix2': 1024.5,
        'crval1': 22.023517632518953,
        'crval2': 11.998755402186378,
        'cd1_1': -1.7436711450380403e-05,
        'cd1_2': -2.0976403747938548e-08,
        'cd2_1': -4.627876291461441e-08,
        'cd2_2': 1.7519986436203685e-05,
        'a_order': 3,
        'b_order': 3,
        'a_0_2': -1.5528350090319145e-06,
        'a_0_3': 4.604566874451819e-13,
        'a_1_1': -1.1606066126743262e-05,
        'a_1_2': 1.712921174539815e-09,
        'a_2_0': 1.9201266376054155e-06,
        'a_2_1': -7.921594730915457e-11,
        'a_3_0': 1.597813188317488e-09,
        'b_0_2': -6.7579614793024975e-06,
        'b_0_3': 1.6914674117189632e-09,
        'b_1_1': 3.633049815338483e-06,
        'b_1_2': -8.645733568868342e-11,
        'b_2_0': 4.914642322124454e-06,
        'b_2_1': 1.5704907984494963e-09,
        'b_3_0': -2.9578595707610732e-12,
        'sipmxerr': 0.1291142984620709,
        'sipiverr': 0.11035960854622291,
        'ap_0_1': 4.186210141016899e-07,
        'ap_0_2': 1.5494761945928307e-06,
        'ap_0_3': 3.783699888145856e-11,
        'ap_1_0': -7.690659272206018e-06,
        'ap_1_1': 1.150139841138833e-05,
        'ap_1_2': -1.499790481486753e-09,
        'ap_2_0': -1.9018787403711958e-06,
        'ap_2_1': -4.510078445255482e-11,
        'ap_3_0': -1.6308694611892058e-09,
        'bp_0_1': -7.464889026127811e-06,
        'bp_0_2': 6.693371868309925e-06,
        'bp_0_3': -1.5900922772608882e-09,
        'bp_1_0': 4.365486723327096e-07,
        'bp_1_1': -3.59966227309358e-06,
        'bp_1_2': -4.470754675266111e-11,
        'bp_2_0': -4.903683300590091e-06,
        'bp_2_1': -1.711464083137801e-09,
        'bp_3_0': 3.894833312782645e-11,
        'ap_order': 3,
        'bp_order': 3
    }

    hdu1 = create_hdul()
    im = ImageModel(hdu1)

    pipe = AssignWcsStep()
    result = pipe.call(im)

    # check that result.meta.wcsinfo has correct
    # values after SIP approx.
    wcs_info = result.meta.wcsinfo.instance

    # make sure all expected keys are there
    assert set(list(true_wcs.keys())).issubset(set((list(wcs_info.keys()))))

    # and make sure they match
    for key in true_wcs:
        if 'ctype' in key:
            assert true_wcs[key] == wcs_info[key]
        else:
            assert_allclose(true_wcs[key], wcs_info[key], rtol=2e-4, atol=1e-10)

    # evaluate fits wcs and gwcs, make sure they agree
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
