import os
import pytest
from astropy.io import fits as pf
from jwst.flatfield.flat_field_step import FlatFieldStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_flat_field_nircam(_bigdata):
    """

    Regression test of flat_field step performed on MIRI data.

    """
    suffix = 'flat_field'
    output_file_base, output_file = add_suffix('flatfield1_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass



    FlatFieldStep.call(_bigdata+'/nircam/test_flat_field/jw00017001001_01101_00001_NRCA1_ramp_fit.fits',
                       output_file=output_file_base, suffix=suffix
                       )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/nircam/test_flat_field/jw00017001001_01101_00001_NRCA1_flat_field.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
