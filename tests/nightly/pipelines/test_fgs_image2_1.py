import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_image2 import Image2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_fgs_image2pipeline1(_bigdata):
    """

    Regression test of calwebb_image2 pipeline performed on FGS imaging mode data.

    """

    Image2Pipeline.call(_bigdata+'/fgs/test_image2pipeline/jw86500007001_02101_00001_GUIDER2_rate.fits',
                        output_file='jw86500007001_02101_00001_GUIDER2_cal.fits')

    h = pf.open('jw86500007001_02101_00001_GUIDER2_cal.fits')
    href = pf.open(_bigdata+'/fgs/test_image2pipeline/jw86500007001_02101_00001_GUIDER2_cal_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['area']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['area']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

