import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_detector1pipeline4(_bigdata):
    """

    Regression test of calwebb_detector1 pipeline performed on NIRSpec data.

    """
    step = Detector1Pipeline()
    step.save_calibrated_ramp = True
    step.ipc.skip = True
    step.persistence.skip = True
    step.jump.rejection_threshold = 4.0
    step.ramp_fit.save_opt = False
    step.output_file = 'jw84600007001_02101_00001_nrs1_rate.fits'
    step.run(_bigdata+'/pipelines/jw84600007001_02101_00001_nrs1_uncal.fits')

    # Compare ramp product
    n_ramp = 'jw84600007001_02101_00001_nrs1_ramp.fits'
    h = fits.open( n_ramp )
    n_ref = _bigdata+'/pipelines/jw84600007001_02101_00001_nrs1_ramp_ref.fits'
    href = fits.open( n_ref )
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],h['pixeldq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare countrate image product
    n_cr = 'jw84600007001_02101_00001_nrs1_rate.fits'
    h = fits.open( n_cr )
    n_ref = _bigdata+'/pipelines/jw84600007001_02101_00001_nrs1_rate_ref.fits'
    href = fits.open( n_ref )
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

