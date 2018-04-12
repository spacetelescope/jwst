import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_detector1pipeline2():
    """

    Regression test of calwebb_detector1 pipeline performed on MIRI data.

    """
    step = Detector1Pipeline()
    step.save_calibrated_ramp = True
    step.refpix.odd_even_columns = True
    step.refpix.use_side_ref_pixels = True
    step.refpix.side_smoothing_length=11
    step.refpix.side_gain=1.0
    step.refpix.odd_even_rows = True
    step.jump.rejection_threshold = 250.0
    step.ramp_fit.save_opt = False

    step.run(BIGDATA+'/miri/test_sloperpipeline/jw80600012001_02101_00003_mirimage_uncal.fits',
             output_file='jw80600012001_02101_00003_mirimage_rate.fits')

    # Compare the calibrated ramp product
    n_cr = 'jw80600012001_02101_00003_mirimage_ramp.fits'
    h = pf.open( n_cr )
    n_ref = BIGDATA+'/miri/test_sloperpipeline/jw80600012001_02101_00003_mirimage_ramp.fits'
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

    # Compare the multi-integration countrate image product
    n_int = 'jw80600012001_02101_00003_mirimage_rateints.fits'
    h = pf.open( n_int )
    n_ref = BIGDATA+'/miri/test_sloperpipeline/jw80600012001_02101_00003_mirimage_rateints.fits'
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

    # Compare the countrate image product
    n_rate = 'jw80600012001_02101_00003_mirimage_rate.fits'
    h = pf.open( n_rate )
    n_ref = BIGDATA+'/miri/test_sloperpipeline/jw80600012001_02101_00003_mirimage_rate.fits'
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

