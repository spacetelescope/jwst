import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_nrs_fs_brightobj_spec2():
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec fixed-slit data
    that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).

    """

    Spec2Pipeline.call(BIGDATA+'/pipelines/jw84600042001_02101_00001_nrs2_rateints.fits',
                       config_file='calwebb_spec2.cfg')

    na = 'jw84600042001_02101_00001_nrs2_calints.fits'
    nb = BIGDATA+'/pipelines/jw84600042001_02101_00001_nrs2_calints_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],h['relsens',1],
                                    h['wavelength',1]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],href['relsens',1],
                                          href['wavelength',1]])
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

    na = 'jw84600042001_02101_00001_nrs2_x1dints.fits'
    nb = BIGDATA+'/pipelines/jw84600042001_02101_00001_nrs2_x1dints_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d',1],h['extract1d',2],h['extract1d',3]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1],href['extract1d',2],href['extract1d',3]])
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

