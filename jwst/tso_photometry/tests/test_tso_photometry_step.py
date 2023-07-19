import pytest

from astropy.io import fits
from datetime import datetime
from jwst.tso_photometry.tso_photometry_step import TSOPhotometryStep


@pytest.fixture
def test_hdu():
    hdu = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['telescop'] = "JWST"
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = datetime.now().date().strftime('%Y-%m-%d')

    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    hdu.append(phdu)
    hdu.append(scihdu)
    return hdu


def test_tsophotometry_process_fails_with_missing_crpix(test_hdu):
    with pytest.raises(ValueError):
        TSOPhotometryStep().process(test_hdu)
