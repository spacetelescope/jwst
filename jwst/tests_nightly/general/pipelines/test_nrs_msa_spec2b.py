import os
import pytest
import shutil

from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_nrs_msa_spec2b(_bigdata):
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec MSA data,
    including barshadow correction.

    """
    file_in = os.path.join(_bigdata, 'pipelines',
                       'jw95065_nrs_msaspec_barshadow.fits')

    step = Spec2Pipeline()
    step.output_file='jw95065_nrs_msaspec_barshadow_cal.fits'
    step.save_bsub = False
    step.save_results = True
    step.resample_spec.skip = True
    step.resample_spec.save_results = True
    step.cube_build.save_results = True
    step.extract_1d.save_results = True
    step.run(file_in)

    # compare _cal file
    na = 'jw95065_nrs_msaspec_barshadow_cal.fits'
    nb = os.path.join(_bigdata,'pipelines',
                      'jw95065_nrs_msaspec_barshadow_cal_ref.fits')
    h = pf.open(na)
    href = pf.open(nb)
    h = h[:-1]
    href = href[:-1]
    result = pf.diff.FITSDiff(h, href,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()

    #na = 'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_s2d.fits'
    #nb = _bigdata+'/pipelines/F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_s2d_ref.fits'
    #h = pf.open(na)
    #href = pf.open(nb)
    #newh = pf.HDUList([h['primary'],h['sci',1],h['wht',1],h['con',1],h['relsens',1],
    #                                h['sci',2],h['wht',2],h['con',2],h['relsens',2],
    #                                h['sci',3],h['wht',3],h['con',3],h['relsens',3],
    #                                h['sci',4],h['wht',4],h['con',4],h['relsens',4],
    #                                h['sci',5],h['wht',5],h['con',5],h['relsens',5]])
    #newhref = pf.HDUList([href['primary'],href['sci',1],href['wht',1],href['con',1],href['relsens',1],
    #                                      href['sci',2],href['wht',2],href['con',2],href['relsens',2],
    #                                      href['sci',3],href['wht',3],href['con',3],href['relsens',3],
    #                                      href['sci',4],href['wht',4],href['con',4],href['relsens',4],
    #                                      href['sci',5],href['wht',5],href['con',5],href['relsens',5]])
    #result = pf.diff.FITSDiff(newh, newhref,
    #                          ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
    #                          rtol = 0.00001)
    assert result.identical, result.report()

    # compare _x1d file
    na = 'jw95065_nrs_msaspec_barshadow_x1d.fits'
    nb = os.path.join(_bigdata,'pipelines',
                      'jw95065_nrs_msaspec_barshadow_x1d_ref.fits')
    h = pf.open(na)
    href = pf.open(nb)
    h = h[:-1]
    href = href[:-1]
    result = pf.diff.FITSDiff(h, href,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()
