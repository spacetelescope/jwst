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


def test_flat_field_nirspec(_bigdata):
    """

    Regression test of flat_field step performed on NIRSpec fixed slit data.

    """
    output_file_base, output_file = add_suffix('flatfield1_output.fits', 'flat_field')

    try:
        os.remove(output_file)
    except:
        pass



    FlatFieldStep.call(_bigdata+'/nirspec/test_flat_field/jw00023001001_01101_00001_NRS1_extract_2d.fits',
                       output_file=output_file_base, name='flat_field'
                       )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/nirspec/test_flat_field/jw00023001001_01101_00001_NRS1_flat_field.fits')
    newh = pf.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],
                                    h['sci',2],h['err',2],h['dq',2],
                                    h['sci',3],h['err',3],h['dq',3],
                                    h['sci',4],h['err',4],h['dq',4]])

    newhref = pf.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],
                                          href['sci',2],href['err',2],href['dq',2],
                                          href['sci',3],href['err',3],href['dq',3],
                                          href['sci',4],href['err',4],href['dq',4]])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
