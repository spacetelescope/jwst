import pytest
from astropy.io import fits 
from jwst.pipeline import Detector1Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_fgs_detector1_1(_bigdata):
    """
    Regression test of calwebb_detector1 pipeline performed on FGS imaging mode data.
    """
    pipe = Detector1Pipeline()
    pipe.ipc.skip = True
    pipe.refpix.odd_even_columns = True
    pipe.refpix.use_side_ref_pixels = True
    pipe.refpix.side_smoothing_length = 11
    pipe.refpix.side_gain = 1.0
    pipe.refpix.odd_even_rows = True
    pipe.jump.rejection_threshold = 250.0
    pipe.persistence.skip = True
    pipe.ramp_fit.save_opt = False
    pipe.save_calibrated_ramp = True
    pipe.output_file = 'jw86500007001_02101_00001_GUIDER2_rate.fits'

    pipe.run(_bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_uncal.fits')

    # Compare calibrated ramp product
    n_cr = 'jw86500007001_02101_00001_GUIDER2_ramp.fits'
    h = fits.open(n_cr)
    n_ref = _bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_ramp_ref.fits'
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['groupdq'],h['pixeldq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['groupdq'],href['pixeldq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare multi-integration countrate image product
    n_int = 'jw86500007001_02101_00001_GUIDER2_rateints.fits'
    h = fits.open(n_int)
    n_ref = _bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_rateints_ref.fits'
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # Compare countrate image product
    n_rate = 'jw86500007001_02101_00001_GUIDER2_rate.fits'
    h = fits.open(n_rate)
    n_ref = _bigdata+'/fgs/test_sloperpipeline/jw86500007001_02101_00001_GUIDER2_rate_ref.fits'
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
