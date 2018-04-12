import os
import pytest
import numpy as np

from numpy.testing import utils

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_miri_image_wcs():
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    output_file_base, output_file = add_suffix('miri_image_wcs_output.fits', 'assignwcsstep')

    try:
        os.remove(output_file)
    except:
        pass

    input_file = os.path.join(BIGDATA, 'miri', 'test_wcs', 'image', 'jw00001001001_01101_00001_MIRIMAGE_ramp_fit.fits')
    ref_file = os.path.join(BIGDATA, 'miri', 'test_wcs', 'image', 'jw00001001001_01101_00001_MIRIMAGE_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file=output_file_base
                       )
    im = ImageModel(output_file)
    imref = ImageModel(ref_file)
    x, y = np.mgrid[:1031, :1024]
    ra, dec = im.meta.wcs(x, y)
    raref, decref = imref.meta.wcs(x, y)
    utils.assert_allclose(ra, raref)
    utils.assert_allclose(dec, decref)
