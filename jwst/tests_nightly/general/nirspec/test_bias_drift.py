import pytest
from astropy.io import fits
from jwst.refpix.refpix_step import RefPixStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_refpix_nirspec(_bigdata):
    """
    Regression test of refpix step performed on NIRSpec data.
    """
    output_file_base, output_file = add_suffix('refpix1_output.fits', 'refpix')

    RefPixStep.call(_bigdata+'/nirspec/test_bias_drift/jw00023001001_01101_00001_NRS1_dq_init.fits',
                    odd_even_columns=True, use_side_ref_pixels=False, side_smoothing_length=10,
                    side_gain=1.0, output_file=output_file_base, name='refpix'
                    )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/nirspec/test_bias_drift/jw00023001001_01101_00001_NRS1_bias_drift.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
