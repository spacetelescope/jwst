import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_image2 import Image2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_image2pipeline1():
    """

    Regression test of calwebb_image2 pipeline performed on MIRI data.

    """

    Image2Pipeline.call(BIGDATA+'/miri/test_image2pipeline/jw00001001001_01101_00001_mirimage_rate.fits'
                        )

    h = pf.open('jw00001001001_01101_00001_mirimage_cal.fits')
    href = pf.open(BIGDATA+'/miri/test_image2pipeline/jw00001001001_01101_00001_MIRIMAGE_uncal_MiriSloperPipeline_Image2Pipeline.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['area']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['area']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

