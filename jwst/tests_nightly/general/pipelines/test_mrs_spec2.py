import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_mrs_spec2():
    """

    Regression test of calwebb_spec2 pipeline performed on MIRI MRS data.

    """

    Spec2Pipeline.call(BIGDATA+'/pipelines/jw10001001001_01101_00001_mirifushort_rate.fits',
                       config_file='calwebb_spec2.cfg', save_results=True)

    na = 'jw10001001001_01101_00001_mirifushort_cal.fits'
    nb = BIGDATA+'/pipelines/jw10001001001_01101_00001_mirifushort_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens2d']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens2d']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jw10001001001_01101_00001_mirifushort_s3d.fits'
    nb = BIGDATA+'/pipelines/jw10001001001_01101_00001_mirifushort_s3d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jw10001001001_01101_00001_mirifushort_x1d.fits'
    nb = BIGDATA+'/pipelines/jw10001001001_01101_00001_mirifushort_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d',1]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
