import os
import pytest
from astropy.io import fits as pf
from jwst.refpix.refpix_step import RefPixStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_refpix_miri2(_bigdata):
    """

    Regression test of refpix step performed on MIRI data.

    """
    suffix = 'refpix'
    output_file_base, output_file = add_suffix('refpix2_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass



    RefPixStep.call(_bigdata+'/miri/test_bias_drift/jw00025001001_01107_00001_MIRIMAGE_saturation.fits',
                    use_side_ref_pixels=False, side_smoothing_length=10, side_gain=1.0,
                    output_file=output_file_base, suffix=suffix
                    )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_bias_drift/jw00025001001_01107_00001_MIRIMAGE_bias_drift.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
