import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_niriss_detector1(_bigdata):
    """
    Regression test of calwebb_detector1 pipeline performed on NIRISS data.
    """
    step = Detector1Pipeline()
    step.save_calibrated_ramp = True
    step.ipc.skip = True
    step.persistence.skip = True
    step.refpix.odd_even_columns = True
    step.refpix.use_side_ref_pixels = True
    step.refpix.side_smoothing_length = 11
    step.refpix.side_gain = 1.0
    step.refpix.odd_even_rows = True
    step.jump.rejection_threshold = 250.0
    step.ramp_fit.save_opt = False
    step.ramp_fit.suffix = 'ramp'
    step.output_file = 'jw00034001001_01101_00001_NIRISS_rate.fits'

    step.run(_bigdata+'/pipelines/jw00034001001_01101_00001_NIRISS_uncal.fits')

    # Compare ramp product
    n_ramp = 'jw00034001001_01101_00001_NIRISS_ramp.fits'
    h = fits.open(n_ramp)
    n_ref = _bigdata+'/pipelines/jw00034001001_01101_00001_NIRISS_ramp_ref.fits'
    href = fits.open( n_ref )
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],href['pixeldq']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    # Compare countrate image product
    n_cr = 'jw00034001001_01101_00001_NIRISS_rate.fits'
    h = fits.open( n_cr )
    n_ref = _bigdata+'/pipelines/jw00034001001_01101_00001_NIRISS_rate_ref.fits'
    href = fits.open( n_ref )
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

