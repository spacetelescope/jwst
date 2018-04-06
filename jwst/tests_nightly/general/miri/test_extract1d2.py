import os
from astropy.io import fits as pf
from jwst.extract_1d.extract_1d_step import Extract1dStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_extract1d_miri2():
    """

    Regression test of extract_1d step performed on MIRI LRS slitless data.

    """
    output_file_base, output_file = add_suffix('extract1d2_output.fits', 'extract_1d')

    try:
        os.remove(output_file)
    except:
        pass



    Extract1dStep.call(BIGDATA+'/miri/test_extract1d/jw80600012001_02101_00003_mirimage_photom.fits',
                       config_file='extract_1d.cfg',
                       output_file=output_file_base
    )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/miri/test_extract1d/jw80600012001_02101_00003_mirimage_x1d.fits')
    newh = pf.HDUList([h['primary'],h['extract1d',1],h['extract1d',2],h['extract1d',3],h['extract1d',4]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1],href['extract1d',2],href['extract1d',3],href['extract1d',4]])
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
