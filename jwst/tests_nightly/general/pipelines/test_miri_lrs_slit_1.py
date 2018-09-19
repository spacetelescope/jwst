import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_miri_lrs_slit_1(_bigdata):
    """
    Regression test of calwebb_spec2 pipeline performed on a single
    MIRI LRS fixed-slit exposure.
    """
    step = Spec2Pipeline()
    step.save_bsub=True,
    step.save_results=True
    step.resample_spec.save_results = True
    step.cube_build.save_results = True
    step.extract_1d.save_results = True
    step.run(_bigdata+'/pipelines/jw00035001001_01101_00001_MIRIMAGE_rate.fits')

    n_cr = 'jw00035001001_01101_00001_MIRIMAGE_cal.fits'
    n_ref = _bigdata+'/pipelines/jw00035001001_01101_00001_MIRIMAGE_cal_ref.fits'
    h = fits.open(n_cr)
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    n_cr = 'jw00035001001_01101_00001_MIRIMAGE_x1d.fits'
    n_ref = _bigdata+'/pipelines/jw00035001001_01101_00001_MIRIMAGE_x1d_ref.fits'
    h = fits.open(n_cr)
    href = fits.open(n_ref)
    newh = fits.HDUList([h['primary'],h['extract1d']])
    newhref = fits.HDUList([href['primary'],href['extract1d']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
