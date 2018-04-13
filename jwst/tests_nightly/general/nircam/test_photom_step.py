import os
import pytest
from astropy.io import fits as pf
from jwst.photom.photom_step import PhotomStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_photom_nircam():
    """

    Regression test of photom step performed on NIRCam imaging data.

    """
    output_file_base, output_file = add_suffix('photom1_output.fits', 'photom')

    try:
        os.remove(output_file)
    except:
        pass



    PhotomStep.call(BIGDATA+'/nircam/test_photom/jw00017001001_01101_00001_NRCA1_emission.fits',
                    output_file=output_file_base
                    )
    h = pf.open(output_file)
    href = pf.open(BIGDATA+'/nircam/test_photom/jw00017001001_01101_00001_NRCA1_photom.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['relsens'],h['area']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['relsens'],href['area']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
