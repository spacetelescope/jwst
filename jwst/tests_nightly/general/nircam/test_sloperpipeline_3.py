import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_detector1pipeline3():
    """

    Regression test of calwebb_detector1 pipeline performed on NIRCam data.

    """

    Detector1Pipeline.call(BIGDATA+'/pipelines/jw82500001003_02101_00001_NRCALONG_uncal.fits',
                        config_file='calwebb_detector1.cfg',
                        output_file='jw82500001003_02101_00001_NRCALONG_rate.fits')

    # Compare ramp product
    n_ramp = 'jw82500001003_02101_00001_NRCALONG_ramp.fits'
    h = pf.open( n_ramp )
    n_ref = BIGDATA+'/pipelines/jw82500001003_02101_00001_NRCALONG_ramp_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],h['pixeldq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the ramp product file - a:', n_ramp )
    print (' ... and the reference file - b:', n_ref)    
   
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare countrate image product
    n_cr = 'jw82500001003_02101_00001_NRCALONG_rate.fits'
    h = pf.open( n_cr )
    n_ref = BIGDATA+'/pipelines/jw82500001003_02101_00001_NRCALONG_rate_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the countrate image product file - a:', n_cr )
    print (' ... and the reference file - b:', n_ref)    
    
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare countrate nints image product
    n_int = 'jw82500001003_02101_00001_NRCALONG_rateints.fits'
    h = pf.open( n_int )
    n_ref = BIGDATA+'/pipelines/jw82500001003_02101_00001_NRCALONG_rateints_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the countrate nints image product file - a:', n_int )
    print (' ... and the reference file - b:', n_ref)
    
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

