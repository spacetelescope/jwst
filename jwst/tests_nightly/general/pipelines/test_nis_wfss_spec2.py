import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_nis_wfss_spec2(_bigdata):
    """
    Regression test of calwebb_spec2 pipeline performed on NIRISS WFSS data.
    """

    Spec2Pipeline.call(_bigdata+'/pipelines/jw87600-a3001_20171109T145456_spec2_001_asn.json')

    na = 'jw87600017001_02101_00002_nis_cal.fits'
    nb = _bigdata+'/pipelines/jw87600017001_02101_00002_nis_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],h['relsens',1],h['wavelength',1],
                                    h['sci',2],h['err',2],h['dq',2],h['relsens',2],
                                    h['sci',3],h['err',3],h['dq',3],h['relsens',3],h['wavelength',3]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],href['relsens',1],href['wavelength',1],
                                          href['sci',2],href['err',2],href['dq',2],href['relsens',2],
                                          href['sci',3],href['err',3],href['dq',3],href['relsens',3],href['wavelength',3]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    na = 'jw87600017001_02101_00002_nis_x1d.fits'
    nb = _bigdata+'/pipelines/jw87600017001_02101_00002_nis_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d',1],h['extract1d',2]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1],href['extract1d',2]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

