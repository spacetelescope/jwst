import pytest
from astropy.io import fits
from jwst.straylight.straylight_step import StraylightStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_straylight1_miri(_bigdata):
    """
    Regression test of straylight performed on MIRI IFUSHORT data.
    """
    suffix = 'straylight'
    output_file_base, output_file = add_suffix('straylight1_output.fits', suffix)

    StraylightStep.call(_bigdata+'/miri/test_straylight/jw80500018001_02101_00002_MIRIFUSHORT_flatfield.fits',
                        output_file=output_file_base, suffix=suffix
                        )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/miri/test_straylight/jw80500018001_02101_00002_MIRIFUSHORT_straylight.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
