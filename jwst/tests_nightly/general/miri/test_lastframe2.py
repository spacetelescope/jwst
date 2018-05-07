import os
import pytest
from astropy.io import fits as pf
from jwst.lastframe.lastframe_step import LastFrameStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_lastframe_miri2(_bigdata):
    """

    Regression test of lastframe step performed on MIRI data.

    """
    suffix = 'lastframe'
    output_file_base, output_file = add_suffix('lastframe2_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass


    LastFrameStep.call(_bigdata+'/miri/test_lastframe/jw80600012001_02101_00003_mirimage_rscd.fits',
                       output_file=output_file_base, suffix=suffix
                       )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_lastframe/jw80600012001_02101_00003_mirimage_lastframe.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
