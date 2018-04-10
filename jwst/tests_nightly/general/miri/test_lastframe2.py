import os
from astropy.io import fits as pf
from jwst.lastframe.lastframe_step import LastFrameStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_lastframe_miri2():
    """

    Regression test of lastframe step performed on MIRI data.

    """
    output_file_base, output_file = add_suffix('lastframe2_output.fits', 'lastframe')

    try:
        os.remove(output_file)
    except:
        pass


    LastFrameStep.call(BIGDATA+'/miri/test_lastframe/jw80600012001_02101_00003_mirimage_rscd.fits',
                       output_file=output_file_base
                       )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/miri/test_lastframe/jw80600012001_02101_00003_mirimage_lastframe.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
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
