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


def test_nirisssoss2pipeline1(_bigdata):
    """

    Regression test of calwebb_spec2 pipeline performed on NIRISS SOSS data.

    """
    collect_pipeline_cfgs()
    args = [
        'calwebb_tso_spec2.cfg',
        op.join(
            _bigdata,
            'pipelines/jw00034001001_01101_00001_NIRISS_rate_ref.fits'
        )
    ]
    Step.from_cmdline(args)

    n_cr = 'jw00034001001_01101_00001_NIRISS_calints.fits'
    h = pf.open(n_cr)
    n_ref = op.join(
        _bigdata,
        'pipelines/jw00034001001_01101_00001_NIRISS_calints_ref.fits'
    )
    href = pf.open(n_ref)
    newh = pf.HDUList([
        h['primary'], h['sci'], h['err'], h['dq'], h['relsens']
    ])
    newhref = pf.HDUList([
        href['primary'], href['sci'], href['err'], href['dq'], href['relsens']
    ])
    result = pf.diff.FITSDiff(
        newh, newhref,
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()

    n_cr = 'jw00034001001_01101_00001_NIRISS_x1dints.fits'
    h = pf.open(n_cr)
    n_ref = _bigdata+'/pipelines/jw00034001001_01101_00001_NIRISS_x1dints_ref.fits'
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'], h['extract1d']])
    newhref = pf.HDUList([href['primary'], href['extract1d']])
    result = pf.diff.FITSDiff(
        newh, newhref,
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()
