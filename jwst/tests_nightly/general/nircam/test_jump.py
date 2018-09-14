import pytest
from astropy.io import fits
from jwst.jump.jump_step import JumpStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_jump_nircam(_bigdata):
    """
    Regression test of jump step performed on NIRCam data.
    """
    suffix = 'jump'
    output_file_base, output_file = add_suffix('jump1_output.fits', suffix)

    JumpStep.call(_bigdata+'/nircam/test_jump/jw00017001001_01101_00001_NRCA1_linearity.fits',
                  rejection_threshold=25.0,
                  output_file=output_file_base, suffix=suffix
                  )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/nircam/test_jump/jw00017001001_01101_00001_NRCA1_jump.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
