import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_nrs_ifu_spec2(_bigdata):
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec IFU data.

    """
    pipe = Spec2Pipeline()
    pipe.save_bsub = True
    pipe.save_results = True
    pipe.resample_spec.save_results = True
    pipe.cube_build.save_results = True
    pipe.extract_1d.save_results = True
    pipe.run(_bigdata+'/pipelines/jw95175001001_02104_00001_nrs1_rate.fits')

    na = 'jw95175001001_02104_00001_nrs1_cal.fits'
    nb = _bigdata+'/pipelines/jw95175001001_02104_00001_nrs1_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens2d'],
                                    h['pathloss_pointsource'],h['wavelength_pointsource'],
                                    h['pathloss_uniformsource'],h['wavelength_uniformsource']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens2d'],
                                    href['pathloss_pointsource'],href['wavelength_pointsource'],
                                    href['pathloss_uniformsource'],href['wavelength_uniformsource']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    na = 'jw95175001001_02104_00001_nrs1_s3d.fits'
    nb = _bigdata+'/pipelines/jw95175001001_02104_00001_nrs1_s3d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    na = 'jw95175001001_02104_00001_nrs1_x1d.fits'
    nb = _bigdata+'/pipelines/jw95175001001_02104_00001_nrs1_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d']])
    newhref = pf.HDUList([href['primary'],href['extract1d']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()
