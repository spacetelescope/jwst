import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_guider import GuiderPipeline


pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_guider_pipeline1(_bigdata):
    """

    Regression test of calwebb_guider pipeline performed on ID-image data.

    """

    GuiderPipeline.call(_bigdata+'/fgs/test_guiderpipeline/jw88600073001_gs-id_7_image-uncal.fits',
                        output_file='jw88600073001_gs-id_7_image-cal.fits')

    # Compare calibrated ramp product
    n_cr = 'jw88600073001_gs-id_7_image-cal.fits'
    h = pf.open( n_cr )
    n_ref = _bigdata+'/fgs/test_guiderpipeline/jw88600073001_gs-id_7_image-cal_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001
    )

    assert result.identical, result.report()


def test_guider_pipeline2(_bigdata):
    """

    Regression test of calwebb_guider pipeline performed on ACQ-1 data.

    """

    GuiderPipeline.call(_bigdata+'/fgs/test_guiderpipeline/jw88600073001_gs-acq1_2016022183837_uncal.fits',
                        output_file='jw88600073001_gs-acq1_2016022183837_cal.fits')

    # Compare calibrated ramp product
    n_cr = 'jw88600073001_gs-acq1_2016022183837_cal.fits'
    h = pf.open( n_cr )
    n_ref = _bigdata+'/fgs/test_guiderpipeline/jw88600073001_gs-acq1_2016022183837_cal_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001
    )

    assert result.identical, result.report()


def test_guider_pipeline3(_bigdata):
    """

    Regression test of calwebb_guider pipeline performed on ID STACKED data.

    """

    GuiderPipeline.call(_bigdata+'/fgs/test_guiderpipeline/jw86600004001_gs-id_1_stacked-uncal.fits',
                        output_file='jw86600004001_gs-id_1_stacked-cal.fits')

    # Compare calibrated ramp product
    n_cr = 'jw86600004001_gs-id_1_stacked-cal.fits'
    h = pf.open(n_cr)
    n_ref = _bigdata+'/fgs/test_guiderpipeline/jw86600004001_gs-id_1_stacked-cal_ref.fits'
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'],h['sci'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001
    )

    assert result.identical, result.report()
