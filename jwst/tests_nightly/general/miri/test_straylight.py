import os
import pytest
from astropy.io import fits as pf
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

    try:
        os.remove(output_file)
    except:
        pass

    StraylightStep.call(_bigdata+'/miri/test_straylight/jw80500018001_02101_00002_MIRIFUSHORT_flatfield.fits',
                        output_file=output_file_base, suffix=suffix
                        )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_straylight/jw80500018001_02101_00002_MIRIFUSHORT_straylight.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
