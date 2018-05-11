import os
import numpy as np
import pytest

from numpy.testing import utils
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_miri_fixed_slit_wcs(_bigdata):
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    output_file_base, output_file = add_suffix('miri_fixed_wcs_output.fits', 'assignwcsstep')

    try:
        os.remove(output_file)
    except:
        pass

    input_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'fixed', 'jw00035001001_01101_00001_mirimage_rate.fits')
    ref_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'fixed', 'jw00035001001_01101_00001_mirimage_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file=output_file_base
                       )
    im = ImageModel(output_file)
    imref = ImageModel(ref_file)
    y, x = np.mgrid[:1031, :1024]
    ra, dec, lam = im.meta.wcs(x, y)
    raref, decref, lamref = imref.meta.wcs(x, y)
    utils.assert_allclose(ra, raref)
    utils.assert_allclose(dec, decref)
    utils.assert_allclose(lam, lamref)
