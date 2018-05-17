import os
import shutil
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_nrs_msa_spec2(_bigdata):
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec MSA data.

    """
    curdir = os.path.abspath(os.curdir)
    # copy over input-specific reference files, such as MSA Metadata file
    msafile = os.path.join(_bigdata, 'pipelines',
                           'jw95065006001_0_short_msa.fits')
    shutil.copy(msafile, curdir)

    input = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'

    # define step for use in test
    step = Spec2Pipeline()
    step.save_bsub = False
    step.output_use_model = True
    step.resample_spec.save_results = True
    step.resample_spec.skip = True
    step.extract_1d.save_results = True
    step.extract_1d.smoothing_length = 0
    step.extract_1d.bkg_order = 0
    step.run(os.path.join(_bigdata, 'pipelines', input))

    output = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod_cal.fits'
    nbname = 'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_cal_ref.fits'
    nb = os.path.join(_bigdata,'pipelines', nbname)
                      
    h = pf.open(output)
    href = pf.open(nb)
    h = h[:-1]
    href = href[:-1]
    result = pf.diff.FITSDiff(h, href,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX','FILENAME'],
                              rtol = 0.00001)

    assert result.identical, result.report()

    output2 = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod_x1d.fits'
    nbname = 'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_x1d_ref.fits'
    nb = os.path.join(_bigdata, 'pipelines', nbname)

    h = pf.open(output2)
    href = pf.open(nb)
    h = h[:-1]
    href = href[:-1]
    result = pf.diff.FITSDiff(h, href,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX','FILENAME'],
                              rtol = 0.00001)

    assert result.identical, result.report()
