import os
import pytest
from astropy.io import fits as pf
from jwst.refpix.refpix_step import RefPixStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_refpix_nirspec_irs2(_bigdata):
    """
    Regression test of refpix step performed on NIRSpec data using IRS2 readout mode.

    """
    output_file_base, output_file = add_suffix('refpix2_output.fits', 'refpix')

    try:
        os.remove(output_file)
    except:
        pass



    RefPixStep.call(_bigdata+'/nirspec/test_bias_drift/jw84600007001_02101_00001_nrs1_superbias.fits',
                    odd_even_columns=True, use_side_ref_pixels=False, side_smoothing_length=10,
                    side_gain=1.0, output_file=output_file_base, name='refpix'
                    )
    h = pf.open(output_file)
    href = pf.open(_bigdata+'/nirspec/test_bias_drift/jw84600007001_02101_00001_nrs1_refpix.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
