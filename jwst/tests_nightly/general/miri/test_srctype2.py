import pytest
from astropy.io import fits
from jwst.srctype.srctype_step import SourceTypeStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_srctype2(_bigdata):
    """
    Regression test of srctype step performed on MIRI LRS slitless data.
    """
    suffix = 'srctyp'
    output_file_base, output_file = add_suffix('srctype2_output.fits', suffix)

    SourceTypeStep.call(_bigdata+'/miri/test_srctype/jw80600012001_02101_00003_mirimage_flat_field.fits',
                        output_file=output_file_base, suffix=suffix
                        )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/miri/test_srctype/jw80600012001_02101_00003_mirimage_srctype.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX', 'FILENAME'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
