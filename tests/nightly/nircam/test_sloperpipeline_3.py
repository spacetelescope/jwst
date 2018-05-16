import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_detector1pipeline3(_bigdata):
    """

    Regression test of calwebb_detector1 pipeline performed on NIRCam data.

    """
    step = Detector1Pipeline()
    step.save_calibrated_ramp = True
    step.ipc.skip = True
    step.refpix.odd_even_columns = True
    step.refpix.use_side_ref_pixels = False
    step.refpix.side_smoothing_length = 10
    step.refpix.side_gain = 1.0
    step.persistence.skip = True
    step.jump.rejection_threshold = 250.0
    step.ramp_fit.save_opt = True
    step.output_file = 'jw82500001003_02101_00001_NRCALONG_rate.fits'
    step.run(_bigdata+'/pipelines/jw82500001003_02101_00001_NRCALONG_uncal.fits')

    # Compare ramp product
    n_ramp = 'jw82500001003_02101_00001_NRCALONG_ramp.fits'
    h = pf.open( n_ramp )
    n_ref = _bigdata+'/pipelines/jw82500001003_02101_00001_NRCALONG_ramp_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],h['pixeldq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare countrate image product
    n_cr = 'jw82500001003_02101_00001_NRCALONG_rate.fits'
    h = pf.open( n_cr )
    n_ref = _bigdata+'/pipelines/jw82500001003_02101_00001_NRCALONG_rate_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare countrate nints image product
    n_int = 'jw82500001003_02101_00001_NRCALONG_rateints.fits'
    h = pf.open( n_int )
    n_ref = _bigdata+'/pipelines/jw82500001003_02101_00001_NRCALONG_rateints_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
