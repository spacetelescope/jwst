import os
from astropy.io import fits as pf
from jwst.superbias.superbias_step import SuperBiasStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_superbias_nirspec():
    """

    Regression test of superbias step performed on NIRSpec data.

    """
    output_file_base, output_file = add_suffix('superbias1_output.fits', 'superbias')

    try:
        os.remove(output_file)
    except:
        pass

    SuperBiasStep.call(BIGDATA+'/nirspec/test_superbias/jw00011001001_01106_00001_NRS2_saturation.fits',
                       config_file='superbias.cfg',
                       output_file=output_file_base
                       )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/nirspec/test_superbias/jw00011001001_01106_00001_NRS2_superbias.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
