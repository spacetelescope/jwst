import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_detector1pipeline1(_bigdata):
    """
    Regression test of calwebb_detector1 pipeline performed on MIRI data.
    """

    step = Detector1Pipeline()
    step.save_calibrated_ramp = True
    step.ipc.skip = True
    step.refpix.odd_even_columns = True
    step.refpix.use_side_ref_pixels = True
    step.refpix.side_smoothing_length=11
    step.refpix.side_gain=1.0
    step.refpix.odd_even_rows = True
    step.persistence.skip = True
    step.jump.rejection_threshold = 250.0
    step.ramp_fit.save_opt = False
    step.output_file='jw00001001001_01101_00001_MIRIMAGE'
    step.suffix='rate'

    step.run(_bigdata+'/miri/test_sloperpipeline/jw00001001001_01101_00001_MIRIMAGE_uncal.fits'
             )

    # Compare calibrated ramp product
    n_cr = 'jw00001001001_01101_00001_MIRIMAGE_ramp.fits'
    h = pf.open( n_cr )
    n_ref = _bigdata+'/miri/test_sloperpipeline/jw00001001001_01101_00001_MIRIMAGE_uncal_jump.fits'
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
    n_int = 'jw00001001001_01101_00001_MIRIMAGE_rateints.fits'
    h = pf.open( n_int )
    n_ref = _bigdata+'/miri/test_sloperpipeline/jw00001001001_01101_00001_MIRIMAGE_uncal_integ.fits'
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
    n_rate = 'jw00001001001_01101_00001_MIRIMAGE_rate.fits'
    h = pf.open( n_rate )
    n_ref = _bigdata+'/miri/test_sloperpipeline/jw00001001001_01101_00001_MIRIMAGE_uncal_MiriSloperPipeline.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
