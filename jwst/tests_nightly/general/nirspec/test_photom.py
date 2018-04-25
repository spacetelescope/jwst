import os
import pytest
from astropy.io import fits as pf
from jwst.photom.photom_step import PhotomStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_photom_nirspec(_bigdata):
    """

    Regression test of photom step performed on NIRSpec fixed slit data.

    """
    output_file_base, output_file = add_suffix('photom1_output.fits', 'photom')

    try:
        os.remove(output_file)
    except:
        pass



    PhotomStep.call(_bigdata+'/nirspec/test_photom/jw00023001001_01101_00001_NRS1_flat_field.fits',
                    output_file=output_file_base, name='photom'
                    )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/nirspec/test_photom/jw00023001001_01101_00001_NRS1_photom.fits')
    newh = pf.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],h['relsens',1],
                                    h['sci',2],h['err',2],h['dq',2],h['relsens',2],
                                    h['sci',3],h['err',3],h['dq',3],h['relsens',3],
                                    h['sci',4],h['err',4],h['dq',4],h['relsens',4]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],href['relsens',1],
                                          href['sci',2],href['err',2],href['dq',2],href['relsens',2],
                                          href['sci',3],href['err',3],href['dq',3],href['relsens',3],
                                          href['sci',4],href['err',4],href['dq',4],href['relsens',4]])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
