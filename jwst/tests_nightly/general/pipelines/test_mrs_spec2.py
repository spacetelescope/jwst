import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_mrs_spec2(_bigdata):
    """

    Regression test of calwebb_spec2 pipeline performed on MIRI MRS data.

    """
    step = Spec2Pipeline()
    step.save_bsub=True,
    step.save_results=True
    step.resample_spec.save_results = True
    step.cube_build.save_results = True
    step.extract_1d.save_results = True
    step.run(_bigdata+'/pipelines/jw10001001001_01101_00001_mirifushort_rate.fits')

    na = 'jw10001001001_01101_00001_mirifushort_cal.fits'
    nb = _bigdata+'/pipelines/jw10001001001_01101_00001_mirifushort_cal_ref.fits'
    h = fits.open(na)
    href = fits.open(nb)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens2d']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens2d']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)
    assert result.identical, result.report()

    na = 'jw10001001001_01101_00001_mirifushort_s3d.fits'
    nb = _bigdata+'/pipelines/jw10001001001_01101_00001_mirifushort_s3d_ref.fits'
    h = fits.open(na)
    href = fits.open(nb)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)
    assert result.identical, result.report()

    na = 'jw10001001001_01101_00001_mirifushort_x1d.fits'
    nb = _bigdata+'/pipelines/jw10001001001_01101_00001_mirifushort_x1d_ref.fits'
    h = fits.open(na)
    href = fits.open(nb)
    newh = fits.HDUList([h['primary'],h['extract1d',1]])
    newhref = fits.HDUList([href['primary'],href['extract1d',1]])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)
    assert result.identical, result.report()
