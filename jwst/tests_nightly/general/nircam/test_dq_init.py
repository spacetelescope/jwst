import os
import pytest
from astropy.io import fits as pf
from jwst.dq_init.dq_init_step import DQInitStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_dq_init_nircam():
    """

    Regression test of dq_init step performed on uncalibrated NIRCam data.

    """
    output_file_base, output_file = add_suffix('dqinit1_output.fits', 'dq_init')

    try:
        os.remove(output_file)
    except:
        pass



    DQInitStep.call(BIGDATA+'/nircam/test_dq_init/jw00017001001_01101_00001_NRCA1_uncal.fits',
                    output_file=output_file_base
                    )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/nircam/test_dq_init/jw00017001001_01101_00001_NRCA1_dq_init.fits')
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
        raise AssertionError(e)
