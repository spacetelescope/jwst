import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

@pytest.mark.xfail(reason='https://github.com/STScI-JWST/jwst/issues/2007')
def test_nis_wfss_spec2(_bigdata):
    """
    Regression test of calwebb_spec2 pipeline performed on NIRISS WFSS data.
    """
    pipe = Spec2Pipeline()
    pipe.save_bsub = True
    pipe.save_results = True
    pipe.resample_spec.save_results = True
    pipe.extract_1d.save_results = True
    pipe.run(_bigdata+'/pipelines/jw87600-a3001_20171109T145456_spec2_001_asn.json')

    # Compare the _cal file
    na = 'jw87600017001_02101_00002_nis_cal.fits'
    nb = _bigdata+'/pipelines/jw87600017001_02101_00002_nis_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    h = h[:-1]
    href = href[:-1]
    result = pf.diff.FITSDiff(h, href,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    # Compare the _x1d file
    na = 'jw87600017001_02101_00002_nis_x1d.fits'
    nb = _bigdata+'/pipelines/jw87600017001_02101_00002_nis_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    h = h[:-1]
    href = href[:-1]
    result = pf.diff.FITSDiff(h, href,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()
