import os
from astropy.io import fits as pf
from jwst.photom.photom_step import PhotomStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_photom_miri2():
    """

    Regression test of photom step performed on MIRI LRS slitless data.

    """
    output_file_base, output_file = add_suffix('photom2_output.fits', 'photom')

    try:
        os.remove(output_file)
    except:
        pass



    PhotomStep.call(BIGDATA+'/miri/test_photom/jw80600012001_02101_00003_mirimage_srctype.fits',
                    config_file='photom.cfg',
                    output_file=output_file_base
    )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/miri/test_photom/jw80600012001_02101_00003_mirimage_photom.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
