import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_image2 import Image2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_image2pipeline2b():
    """

    Regression test of calwebb_image2 pipeline performed on NIRCam data,
    using a multiple integration rate (rateints) file as input.

    """

    Image2Pipeline.call(BIGDATA+'/pipelines/jw82500001003_02101_00001_NRCALONG_rateints.fits',
                        output_file='jw82500001003_02101_00001_NRCALONG_calints.fits')

    h = pf.open('jw82500001003_02101_00001_NRCALONG_calints.fits')
    href = pf.open(BIGDATA+'/pipelines/jw82500001003_02101_00001_NRCALONG_calints_ref.fits')
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

