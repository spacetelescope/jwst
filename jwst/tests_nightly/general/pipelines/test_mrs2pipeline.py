import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_mrs2pipeline1():
    """

    Regression test of calwebb_spec2 pipeline performed on MIRI MRS data.

    """

    Spec2Pipeline.call(BIGDATA+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_rate.fits',
                       config_file='calwebb_spec2.cfg', save_results=True)

    n_h = 'jw80500018001_02101_00002_MIRIFUSHORT_cal.fits'
    h = pf.open(n_h)
    n_href = BIGDATA+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_cal.fits'
    href = pf.open(n_href)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between the calibrated product file - a:', n_h)
    print (' ... and the reference file - b:', n_href)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    n_h = 'jw80500018001_02101_00002_MIRIFUSHORT_s3d.fits'
    h = pf.open(n_h)
    n_href = BIGDATA+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_s3d.fits'
    href = pf.open(n_href)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between the IFU cube product file - a:', n_h)
    print (' ... and the reference file - b:', n_href)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    n_h = 'jw80500018001_02101_00002_MIRIFUSHORT_x1d.fits'
    h = pf.open(n_h)
    n_href = BIGDATA+'/miri/test_mrs2pipeline/jw80500018001_02101_00002_MIRIFUSHORT_x1d.fits'
    href = pf.open(n_href)
    newh = pf.HDUList([h['primary'],h['extract1d']])
    newhref = pf.HDUList([href['primary'],href['extract1d']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between the extracted 1d product file - a:', n_h)
    print (' ... and the reference file - b:', n_href)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

