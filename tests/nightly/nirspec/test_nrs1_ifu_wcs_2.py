import os
import pytest
import numpy as np
from numpy.testing import utils

from jwst.assign_wcs import AssignWcsStep, nirspec

from jwst.datamodels import ImageModel

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_nirspec_nrs1_wcs(_bigdata):
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    output_file_base, output_file = add_suffix('nrs1_ifu_wcs_output.fits', 'assignwcsstep')

    try:
        os.remove(output_file)
    except:
        pass

    input_file = os.path.join(_bigdata, 'nirspec', 'test_wcs', 'nrs1-ifu', 'jw00011001001_01120_00001_NRS1_rate_opaque.fits')
    ref_file = os.path.join(_bigdata, 'nirspec', 'test_wcs', 'nrs1-ifu', 'jw00011001001_01120_00001_NRS1_rate_opaque_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file=output_file_base, name='assignwcsstep'
                       )
    im = ImageModel(output_file)
    imref = ImageModel(ref_file)
    a_wcs = nirspec.nrs_wcs_set_input(im, 0)
    w = a_wcs
    y, x = np.mgrid[w.bounding_box[1][0]:w.bounding_box[1][1], w.bounding_box[0][0]: w.bounding_box[0][1]]
    ra, dec, lam = w(x, y)
    a_wcs_ref = nirspec.nrs_wcs_set_input(im, 0)
    wref = a_wcs_ref
    raref, decref, lamref = wref(x, y)

    # equal_nan is used here as many of the entries are nan.
    # The domain is defined but it is only a few entries in there that are valid
    # as it is a curved narrow slit.
    utils.assert_allclose(ra, raref, equal_nan=True)
    utils.assert_allclose(dec, decref, equal_nan=True)
    utils.assert_allclose(lam, lamref, equal_nan=True)
