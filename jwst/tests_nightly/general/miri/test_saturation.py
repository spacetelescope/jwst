import os
import pytest
from astropy.io import fits as pf
from jwst.saturation.saturation_step import SaturationStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_saturation_miri(_bigdata):
    """

    Regression test of saturation step performed on uncalibrated MIRI data.

    """
    suffix = 'saturation'
    output_file_base, output_file = add_suffix('saturation1_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass



    SaturationStep.call(_bigdata+'/miri/test_saturation/jw00001001001_01101_00001_MIRIMAGE_dq_init.fits',
                        output_file=output_file_base, suffix=suffix
                        )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_saturation/jw00001001001_01101_00001_MIRIMAGE_saturation.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
