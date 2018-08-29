import os
import pytest
from astropy.io import fits as pf
from jwst.wfs_combine.wfs_combine_step import WfsCombineStep

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_wfs_combine(_bigdata):
    """

    Regression test of wfs_combine using do_refine=False (default)
    Association table has 3 (identical) pairs of input files to combine

    """

    WfsCombineStep.call(_bigdata+'/nircam/test_wfs_combine/wfs_3sets_asn.json')

    # compare 1st pair of output files
    h = pf.open('test_wfscom_wfscmb.fits')
    href = pf.open(_bigdata+'/nircam/test_wfs_combine/test_wfscom.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()


    # compare 2nd pair of output files
    h = pf.open('test_wfscoma_wfscmb.fits')
    href = pf.open(_bigdata+'/nircam/test_wfs_combine/test_wfscoma.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # compare 3rd pair of output files
    h = pf.open('test_wfscomb_wfscmb.fits')
    href = pf.open(_bigdata+'/nircam/test_wfs_combine/test_wfscomb.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
