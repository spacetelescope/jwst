import pytest
from astropy.io import fits
from jwst.dark_current.dark_current_step import DarkCurrentStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_dark_current_miri(_bigdata):
    """
    Regression test of dark current step performed on MIRI data.
    """
    suffix = 'dark_current'
    output_file_base, output_file = add_suffix('darkcurrent1_output.fits', suffix)

    DarkCurrentStep.call(_bigdata+'/miri/test_dark_step/jw00001001001_01101_00001_MIRIMAGE_bias_drift.fits',
                         output_file=output_file_base, suffix=suffix
                         )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/miri/test_dark_step/jw00001001001_01101_00001_MIRIMAGE_dark_current.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
