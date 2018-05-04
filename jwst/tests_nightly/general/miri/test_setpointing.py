import os
from shutil import copyfile
from astropy.io import fits as pf
from jwst.lib.set_telescope_pointing import add_wcs

BIGDATA = os.environ['TEST_BIGDATA']

def test_miri_setpointing():
    """
    Regression test of the set_telescope_pointing script on a level-1b MIRI file.
    """

    try:
        os.remove(BIGDATA+'/miri/test_pointing/jw80600010001_02101_00001_mirimage_uncal.fits')
    except:
        pass

    # Copy original version of file to test file, which will get overwritten by test
    copyfile(BIGDATA+'/miri/test_pointing/jw80600010001_02101_00001_mirimage_uncal_orig.fits',
             BIGDATA+'/miri/test_pointing/jw80600010001_02101_00001_mirimage_uncal.fits')

    add_wcs(BIGDATA+'/miri/test_pointing/jw80600010001_02101_00001_mirimage_uncal.fits')

    h = pf.open(BIGDATA+'/miri/test_pointing/jw80600010001_02101_00001_mirimage_uncal.fits')
    href = pf.open(BIGDATA+'/miri/test_pointing/jw80600010001_02101_00001_mirimage_uncal_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['refout'],h['group']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['refout'],href['group']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
