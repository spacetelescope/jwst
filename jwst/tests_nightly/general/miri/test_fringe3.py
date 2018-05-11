import os
import pytest
from astropy.io import fits as pf
from jwst.fringe.fringe_step import FringeStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_fringe_miri3(_bigdata):
    """

    Regression test of fringe performed on MIRI data.

    """
    suffix = 'fringe'
    output_file_base, output_file = add_suffix('fringe3_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass

    FringeStep.call(_bigdata+'/miri/test_fringe/fringe3_input.fits',
                    output_file=output_file_base, suffix=suffix
                    )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_fringe/baseline_fringe3.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
