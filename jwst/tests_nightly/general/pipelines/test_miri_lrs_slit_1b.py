import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_miri_lrs_slit_1b():
    """

    Regression test of calwebb_spec2 pipeline performed on a single
    MIRI LRS fixed-slit exposure with multiple integrations.

    """

    Spec2Pipeline.call(BIGDATA+'/pipelines/jw00035001001_01101_00001_MIRIMAGE_rateints.fits',
                       config_file='calwebb_spec2.cfg')

    n_cr = 'jw00035001001_01101_00001_MIRIMAGE_calints.fits'
    n_ref = BIGDATA+'/pipelines/jw00035001001_01101_00001_MIRIMAGE_calints_ref.fits'
    h = pf.open(n_cr)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the product file - a:', n_cr )
    print (' ... and the reference file - b:', n_ref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    n_cr = 'jw00035001001_01101_00001_MIRIMAGE_x1dints.fits'
    n_ref = BIGDATA+'/pipelines/jw00035001001_01101_00001_MIRIMAGE_x1dints_ref.fits'
    h = pf.open(n_cr)
    href = pf.open(n_ref)
    newh = pf.HDUList([h['primary'],h['extract1d',1],h['extract1d',2],h['extract1d',3],h['extract1d',4]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1],href['extract1d',2],href['extract1d',3],href['extract1d',4]])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the product file - a:', n_cr )
    print (' ... and the reference file - b:', n_ref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

