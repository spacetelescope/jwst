import os
import pytest
from astropy.io import fits as pf
from jwst.superbias.superbias_step import SuperBiasStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_superbias_nirspec(_bigdata):
    """

    Regression test of superbias step performed on NIRSpec data.

    """
    output_file_base, output_file = add_suffix('superbias1_output.fits', 'superbias')

    try:
        os.remove(output_file)
    except:
        pass

    SuperBiasStep.call(_bigdata+'/nirspec/test_superbias/jw00011001001_01106_00001_NRS2_saturation.fits',
                       output_file=output_file_base, name='superbias'
                       )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/nirspec/test_superbias/jw00011001001_01106_00001_NRS2_superbias.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
