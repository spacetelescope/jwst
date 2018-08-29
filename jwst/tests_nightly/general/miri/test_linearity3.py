import os
import pytest
from astropy.io import fits as pf
from jwst.linearity.linearity_step import LinearityStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_linearity_miri3(_bigdata):
    """

    Regression test of linearity step performed on MIRI data.

    """
    suffix = 'linearity'
    output_file_base, output_file = add_suffix('linearity3_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass



    LinearityStep.call(_bigdata+'/miri/test_linearity/jw00001001001_01109_00001_MIRIMAGE_dark_current.fits',
                       override_linearity=_bigdata+'/miri/test_linearity/lin_nan_flag_miri.fits',
                       output_file=output_file_base, suffix=suffix
                       )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_linearity/jw00001001001_01109_00001_MIRIMAGE_linearity.fits')

    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
