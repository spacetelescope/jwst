import pytest
from astropy.io import fits
from jwst.photom.photom_step import PhotomStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_photom_miri(_bigdata):
    """
    Regression test of photom step performed on MIRI imaging data.
    """
    suffix = 'photom'
    output_file_base, output_file = add_suffix('photom1_output.fits', suffix)

    PhotomStep.call(_bigdata+'/miri/test_photom/jw00001001001_01101_00001_MIRIMAGE_emission.fits',
                    output_file=output_file_base, suffix=suffix
                    )
    h = fits.open(output_file)
    href = fits.open(_bigdata+'/miri/test_photom/jw00001001001_01101_00001_MIRIMAGE_photom.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens'],h['area']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens'],href['area']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
