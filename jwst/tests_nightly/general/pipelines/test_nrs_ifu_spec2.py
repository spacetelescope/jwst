import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_nrs_ifu_spec2():
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec IFU data.

    """

    Spec2Pipeline.call(BIGDATA+'/pipelines/jw95175001001_02104_00001_nrs1_rate.fits',
                       config_file='calwebb_spec2.cfg')

    na = 'jw95175001001_02104_00001_nrs1_cal.fits'
    nb = BIGDATA+'/pipelines/jw95175001001_02104_00001_nrs1_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens2d'],
                                    h['pathloss_pointsource'],h['wavelength_pointsource'],
                                    h['pathloss_uniformsource'],h['wavelength_uniformsource']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens2d'],
                                    href['pathloss_pointsource'],href['wavelength_pointsource'],
                                    href['pathloss_uniformsource'],href['wavelength_uniformsource']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb) 

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jw95175001001_02104_00001_nrs1_s3d.fits'
    nb = BIGDATA+'/pipelines/jw95175001001_02104_00001_nrs1_s3d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb) 

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jw95175001001_02104_00001_nrs1_x1d.fits'
    nb = BIGDATA+'/pipelines/jw95175001001_02104_00001_nrs1_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d']])
    newhref = pf.HDUList([href['primary'],href['extract1d']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb) 

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

