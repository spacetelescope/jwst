import pytest
from astropy.io import fits
from jwst.extract_1d.extract_1d_step import Extract1dStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_extract1d_miri(_bigdata):
    """
    Regression test of extract_1d step performed on MIRI LRS fixed-slit data.
    """
    suffix = 'extract_1d'
    output_file_base, output_file = add_suffix('extract1d1_output.fits', suffix)

    Extract1dStep.call(_bigdata+'/miri/test_extract1d/jw00035001001_01101_00001_mirimage_photom.fits',
                       smoothing_length=0,
                       output_file=output_file_base, suffix=suffix
                       )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/miri/test_extract1d/jw00035001001_01101_00001_mirimage_x1d.fits')
    newh = fits.HDUList([h['primary'],h['extract1d',1]])
    newhref = fits.HDUList([href['primary'],href['extract1d',1]])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
