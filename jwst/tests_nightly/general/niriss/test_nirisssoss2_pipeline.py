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


def test_nirisssoss2pipeline1(_bigdata):
    """
    Regression test of calwebb_tso-spec2 pipeline performed on NIRISS SOSS data.
    """
    collect_pipeline_cfgs()
    args = [
        'calwebb_tso-spec2.cfg',
        op.join(
            _bigdata,
            'pipelines/jw10003001002_03101_00001-seg003_nis_rateints.fits'
        )
    ]
    Step.from_cmdline(args)

    # Compare the _calints products
    n_cr = 'jw10003001002_03101_00001-seg003_nis_calints.fits'
    n_ref = _bigdata+'/pipelines/jw10003001002_03101_00001-seg003_nis_calints_ref.fits'
    h = fits.open(n_cr)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(
        h, href,
        ignore_hdus=['INT_TIMES', 'VAR_POISSON', 'VAR_RNOISE', 'ASDF'],
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()

    # Compare the _x1dints products
    n_cr = 'jw10003001002_03101_00001-seg003_nis_x1dints.fits'
    n_ref = _bigdata+'/pipelines/jw10003001002_03101_00001-seg003_nis_x1dints_ref.fits'
    h = fits.open(n_cr)
    href = fits.open(n_ref)
    result = fits.diff.FITSDiff(
        h, href,
        ignore_hdus=['INT_TIMES', 'ASDF'],
        ignore_keywords=['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX'],
        rtol=0.00001
    )
    assert result.identical, result.report()
