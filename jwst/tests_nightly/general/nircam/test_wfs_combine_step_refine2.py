import os
import pytest
from astropy.io import fits as pf
from jwst.wfs_combine.wfs_combine_step import WfsCombineStep 

BIGDATA = os.environ['TEST_BIGDATA']

def test_wfs_combine():
    """

    Regression test of wfs_combine using do_refine=True 

    """

    WfsCombineStep.call(BIGDATA+'/nircam/test_wfs_combine/wfs_3sets_asn3.json',
                         config_file='wfs_combine2.cfg' )

    # compare 1st pair of output files
    h = pf.open('test_wfscom3_wfscmb.fits')
    href = pf.open(BIGDATA+'/nircam/test_wfs_combine/test_wfscom_do_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX','FILENAME'],
                              rtol = 0.00001
    )
    result.report()
 
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # compare 2nd pair of output files
    h = pf.open('test_wfscom3a_wfscmb.fits')
    href = pf.open(BIGDATA+'/nircam/test_wfs_combine/test_wfscoma_do_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX','FILENAME'],
                              rtol = 0.00001
    )
    result.report()
 
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # compare 3rd pair of output files
    h = pf.open('test_wfscom3b_wfscmb.fits')
    href = pf.open(BIGDATA+'/nircam/test_wfs_combine/test_wfscomb_do_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX','FILENAME'],
                              rtol = 0.00001
    )
    result.report()
 
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
