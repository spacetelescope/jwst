import os
import pytest
from astropy.io import fits
from jwst.extract_1d.extract_1d_step import Extract1dStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_extract1d_nirspec(_bigdata):
    """Regression test of extract_1d step performed on NIRSpec fixed slit data.
    """
    output_file_base, output_file = add_suffix('extract1d1_output.fits', 'extract_1d')

    try:
        os.remove(output_file)
    except:
        pass

    Extract1dStep.call(_bigdata + '/nirspec/test_extract_1d/jw00023001001_01101_00001_NRS1_cal.fits',
                       smoothing_length=0, bkg_order=0,
                       output_file=output_file_base, name='extract_1d')

    h = fits.open(output_file)
    href = fits.open(_bigdata + '/nirspec/test_extract_1d/jw00023001001_01101_00001_NRS1_spec.fits')
    newh = fits.HDUList([h['primary'],
                         h[('extract1d',1)],
                         h[('extract1d',2)],
                         h[('extract1d',3)],
                         h[('extract1d',4)]])
    newhref = fits.HDUList([href['primary'],
                            href[('extract1d',1)],
                            href[('extract1d',2)],
                            href[('extract1d',3)],
                            href[('extract1d',4)]])
    result = fits.diff.FITSDiff(newh,
                                newhref,
                                ignore_keywords=['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                                rtol=0.00001)
    assert result.identical, result.report()
