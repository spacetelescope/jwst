import os
import pytest

import numpy as np
from numpy.testing import utils

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import CubeModel

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_miri_slitless_wcs(_bigdata):
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    suffix = 'assignwcsstep'
    output_file_base, output_file = add_suffix('miri_slitless_wcs_output.fits', suffix)

    try:
        os.remove(output_file)
    except:
        pass

    input_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'slitless', 'jw80600012001_02101_00003_mirimage_rateints.fits')
    ref_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'slitless', 'jw80600012001_02101_00003_mirimage_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file=output_file_base, suffix=suffix
                       )
    im = CubeModel(output_file)
    imref = CubeModel(ref_file)
    x, y = np.mgrid[:1031, :1024]
    ra, dec, lam = im.meta.wcs(x, y)
    raref, decref, lamref = imref.meta.wcs(x, y)
    utils.assert_allclose(ra, raref)
    utils.assert_allclose(dec, decref)
    utils.assert_allclose(lam, lamref)
