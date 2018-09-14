from os import path as op
import pytest
from astropy.io import fits

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_mirilrs2pipeline1(_bigdata):
    """
    Regression test of calwebb_spec2 pipeline performed on
    MIRI LRS slitless data.
    """
    collect_pipeline_cfgs()
    args = [
        'calwebb_tso_spec2.cfg',
        op.join(
            _bigdata,
            'pipelines/jw80600012001_02101_00003_mirimage_rateints.fits'
        ),
    ]
    Step.from_cmdline(args)

    n_cr = 'jw80600012001_02101_00003_mirimage_calints.fits'
    h = fits.open(n_cr)
    n_ref = op.join(
        _bigdata,
        'pipelines/jw80600012001_02101_00003_mirimage_calints_ref.fits'
    )
    href = fits.open(n_ref)
    newh = fits.HDUList([
        h['primary'], h['sci'], h['err'], h['dq'], h['relsens']
    ])
    newhref = fits.HDUList([
        href['primary'], href['sci'], href['err'], href['dq'], href['relsens']
    ])
    result = fits.diff.FITSDiff(
        newh, newhref,
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()

    n_cr = 'jw80600012001_02101_00003_mirimage_x1dints.fits'
    h = fits.open(n_cr)
    n_ref = op.join(
        _bigdata,
        'pipelines/jw80600012001_02101_00003_mirimage_x1dints_ref.fits'
    )
    href = fits.open(n_ref)
    newh = fits.HDUList([
        h['primary'], h['extract1d', 1], h['extract1d', 2], h['extract1d', 3], h['extract1d', 4]
    ])
    newhref = fits.HDUList([
        href['primary'], href['extract1d', 1], href['extract1d', 2], href['extract1d', 3], href['extract1d', 4]
    ])
    result = fits.diff.FITSDiff(
        newh, newhref,
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()
