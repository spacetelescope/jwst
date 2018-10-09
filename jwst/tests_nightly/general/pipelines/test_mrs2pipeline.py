import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_mrs2pipeline1(_bigdata):
    """

    Regression test of calwebb_spec2 pipeline performed on MIRI MRS data.

    """
    step = Spec2Pipeline()
    step.save_bsub=True,
    step.save_results=True
    step.resample_spec.save_results = True
    step.cube_build.save_results = True
    step.extract_1d.save_results = True
    step.run(_bigdata+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_rate.fits')

    n_h = 'jw80500018001_02101_00002_MIRIFUSHORT_cal.fits'
    h = fits.open(n_h)
    n_href = _bigdata+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_cal.fits'
    href = fits.open(n_href)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    n_h = 'jw80500018001_02101_00002_MIRIFUSHORT_s3d.fits'
    h = fits.open(n_h)
    n_href = _bigdata+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_s3d.fits'
    href = fits.open(n_href)
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    n_h = 'jw80500018001_02101_00002_MIRIFUSHORT_x1d.fits'
    h = fits.open(n_h)
    n_href = _bigdata+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_x1d.fits'
    href = fits.open(n_href)
    newh = fits.HDUList([h['primary'],h['extract1d']])
    newhref = fits.HDUList([href['primary'],href['extract1d']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()
