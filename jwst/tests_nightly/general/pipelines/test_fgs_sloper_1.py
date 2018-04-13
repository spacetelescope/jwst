import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_fgs_detector1_1():
    """

    Regression test of calwebb_detector1 pipeline performed on FGS imaging mode data.

    """

    Detector1Pipeline.call(_bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_uncal.fits',
                           save_calibrated_ramp=True,
                           output_file='jw86500007001_02101_00001_GUIDER2_rate.fits')

    # Compare calibrated ramp product
    n_cr = 'jw86500007001_02101_00001_GUIDER2_ramp.fits'
    h = pf.open( n_cr )
    n_ref = _bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_ramp_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],href['pixeldq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare multi-integration countrate image product
    n_int = 'jw86500007001_02101_00001_GUIDER2_rateints.fits'
    h = pf.open( n_int )
    n_ref = _bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_rateints_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare countrate image product
    n_rate = 'jw86500007001_02101_00001_GUIDER2_rate.fits'
    h = pf.open( n_rate )
    n_ref = _bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_rate_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
