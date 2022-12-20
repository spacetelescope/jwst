"""
Test script for set_velocity_aberration.py
"""
import subprocess

from numpy import isclose
from astropy.io import fits

from jwst.lib.set_velocity_aberration import compute_va_effects

# Testing constants
GOOD_VELOCITY = (100.0, 100.0, 100.0)
GOOD_POS = (359.0, -2.0)
GOOD_SCALE_FACTOR = 1.000316017905845
GOOD_APPARENT_RA = 359.01945099823
GOOD_APPARENT_DEC = -1.980247580394956


def test_compute_va_effects_valid():
    scale_factor, va_ra, va_dec = compute_va_effects(*GOOD_VELOCITY, *GOOD_POS)
    assert isclose(scale_factor, GOOD_SCALE_FACTOR)
    assert isclose(va_ra, GOOD_APPARENT_RA)
    assert isclose(va_dec, GOOD_APPARENT_DEC)


def test_compute_va_effects_zero_velocity():
    scale_factor, va_ra, va_dec = compute_va_effects(0.0, 0.0, 0.0, *GOOD_POS)
    assert isclose(scale_factor, 1.0, atol=1e-16)
    assert isclose(va_ra, GOOD_POS[0], atol=1e-16)
    assert isclose(va_dec, GOOD_POS[1], atol=1e-16)


def test_velocity_aberration_script(tmpdir):
    """Test the whole script on a FITS file"""
    path = str(tmpdir.join("velocity_aberration_tmpfile.fits"))
    hdulist = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(name='SCI')])
    hdulist['SCI'].header['JWST_DX'] = GOOD_VELOCITY[0]
    hdulist['SCI'].header['JWST_DY'] = GOOD_VELOCITY[1]
    hdulist['SCI'].header['JWST_DZ'] = GOOD_VELOCITY[2]
    hdulist['SCI'].header['RA_REF'] = GOOD_POS[0]
    hdulist['SCI'].header['DEC_REF'] = GOOD_POS[1]
    hdulist.writeto(path)
    hdulist.close()

    subprocess.check_call(["set_velocity_aberration", path])

    with fits.open(path) as hdulist_in:
        assert isclose(hdulist_in[0].header['VA_RA'], GOOD_APPARENT_RA, rtol=0, atol=1e-7)
        assert isclose(hdulist_in[0].header['VA_DEC'], GOOD_APPARENT_DEC, rtol=0, atol=1e-7)
        assert isclose(hdulist_in['SCI'].header['VA_SCALE'], GOOD_SCALE_FACTOR, rtol=0, atol=1e-7)
