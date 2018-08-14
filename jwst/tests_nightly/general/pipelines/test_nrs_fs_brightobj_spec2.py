from os import path as op
import pytest
from astropy.io import fits as pf

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_nrs_fs_brightobj_spec2(_bigdata):
    """
    Regression test of calwebb_spec2 pipeline performed on NIRSpec
    fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).
    """
    collect_pipeline_cfgs()
    args = [
        'calwebb_tso_spec2.cfg',
        op.join(
            _bigdata,
            'pipelines/jw84600042001_02101_00001_nrs2_rateints.fits'
        )
    ]
    Step.from_cmdline(args)

    ignore_keywords = [
        'DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'
    ]

    na = 'jw84600042001_02101_00001_nrs2_calints.fits'
    nb = op.join(
        _bigdata,
        'pipelines/jw84600042001_02101_00001_nrs2_calints_ref.fits'
    )
    h = pf.open(na)
    href = pf.open(nb)
    result = pf.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.00001)
    assert result.identical, result.report()

    na = 'jw84600042001_02101_00001_nrs2_x1dints.fits'
    nb = op.join(
        _bigdata,
        'pipelines/jw84600042001_02101_00001_nrs2_x1dints_ref.fits'
    )
    h = pf.open(na)
    href = pf.open(nb)
    result = pf.diff.FITSDiff(h,
                              href,
                              ignore_hdus=['ASDF'],
                              ignore_keywords=ignore_keywords,
                              rtol=0.00001)
    assert result.identical, result.report()
