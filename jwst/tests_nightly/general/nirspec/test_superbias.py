import pytest
from astropy.io import fits

from jwst.superbias import SuperBiasStep

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

    SuperBiasStep.call(_bigdata+'/nirspec/test_superbias/jw00011001001_01106_00001_NRS2_saturation.fits',
                       output_file=output_file_base, name='superbias'
                       )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/nirspec/test_superbias/jw00011001001_01106_00001_NRS2_superbias.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
