import os
from astropy.io import fits as pf
from jwst.fringe.fringe_step import FringeStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_fringe_miri3():
    """

    Regression test of fringe performed on MIRI data.

    """
    output_file_base, output_file = add_suffix('fringe3_output.fits', 'fringe')

    try:
        os.remove(output_file)
    except:
        pass

    FringeStep.call(BIGDATA+'/miri/test_fringe/fringe3_input.fits',
                    output_file=output_file_base
                    )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/miri/test_fringe/baseline_fringe3.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
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
