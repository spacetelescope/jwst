import os
import pytest
from astropy.io import fits as pf
from jwst.photom.photom_step import PhotomStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_photom_miri2(_bigdata):
    """

    Regression test of photom step performed on MIRI LRS slitless data.

    """
    suffix = 'photom'
    output_file_base, output_file = add_suffix('photom2_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass



    PhotomStep.call(_bigdata+'/miri/test_photom/jw80600012001_02101_00003_mirimage_srctype.fits',
                    output_file=output_file_base, suffix=suffix
                    )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_photom/jw80600012001_02101_00003_mirimage_photom.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
