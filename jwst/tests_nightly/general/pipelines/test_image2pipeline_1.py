import pytest

from astropy.io import fits

from jwst.pipeline.calwebb_image2 import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_image2pipeline1(_bigdata):
    """
    Regression test of calwebb_image2 pipeline performed on MIRI data.
    """

    collect_pipeline_cfgs('cfgs')
    Image2Pipeline.call(_bigdata+'/miri/test_image2pipeline/jw00001001001_01101_00001_mirimage_rate.fits',
                        config_file='cfgs/calwebb_image2.cfg'
    )

    h = fits.open('jw00001001001_01101_00001_mirimage_cal.fits')
    href = fits.open(_bigdata+'/miri/test_image2pipeline/jw00001001001_01101_00001_mirimage_cal_ref.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['area']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['area']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    h = fits.open('jw00001001001_01101_00001_mirimage_i2d.fits')
    href = fits.open(_bigdata+'/miri/test_image2pipeline/jw00001001001_01101_00001_mirimage_i2d_ref.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['wht'],h['con']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['wht'],href['con']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
