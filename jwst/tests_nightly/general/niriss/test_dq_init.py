import pytest
from astropy.io import fits
from jwst.dq_init.dq_init_step import DQInitStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_dq_init_niriss(_bigdata):
    """
    Regression test of dq_init step performed on uncalibrated NIRISS data.
    """
    suffix = 'dq_init'
    output_file_base, output_file = add_suffix('dqinit1_output.fits', suffix)

    DQInitStep.call(_bigdata+'/niriss/test_dq_init/jw00034001001_01101_00001_NIRISS_uncal.fits',
                    output_file=output_file_base, suffix=suffix
                    )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/niriss/test_dq_init/jw00034001001_01101_00001_NIRISS_dq_init.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
