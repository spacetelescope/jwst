import pytest
from astropy.io import fits
from jwst.emission.emission_step import EmissionStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_emission_nircam(_bigdata):
    """
    Regression test of emission step performed on calibrated NIRCam data.
    """
    suffix = 'emission'
    output_file_base, output_file = add_suffix('emission1_output.fits', suffix)

    EmissionStep.call(_bigdata+'/nircam/test_emission/jw00017001001_01101_00001_NRCA1_persistence.fits',
                      output_file=output_file_base, suffix=suffix
                      )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/nircam/test_emission/jw00017001001_01101_00001_NRCA1_emission.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
