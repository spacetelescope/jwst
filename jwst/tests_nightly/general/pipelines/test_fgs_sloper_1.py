import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_fgs_detector1_1():
    """

    Regression test of calwebb_detector1 pipeline performed on FGS imaging mode data.

    """

    Detector1Pipeline.call(BIGDATA+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_uncal.fits',
                           save_calibrated_ramp=True,
                           output_file='jw86500007001_02101_00001_GUIDER2_rate.fits')

    # Compare calibrated ramp product
    n_cr = 'jw86500007001_02101_00001_GUIDER2_ramp.fits'
    h = pf.open( n_cr )
    n_ref = BIGDATA+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_ramp_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],href['pixeldq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the calibrated ramp product file - a:', n_cr )
    print (' ... and the reference file - b:', n_ref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare multi-integration countrate image product
    n_int = 'jw86500007001_02101_00001_GUIDER2_rateints.fits'
    h = pf.open( n_int )
    n_ref = BIGDATA+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_rateints_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the multi-integration countrate image product file - a:', n_int )
    print (' ... and the reference file - b:', n_ref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # Compare countrate image product
    n_rate = 'jw86500007001_02101_00001_GUIDER2_rate.fits'
    h = pf.open( n_rate )
    n_ref = BIGDATA+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_rate_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the countrate image product file - a:', n_rate )
    print (' ... and the reference file - b:', n_ref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

