import os
import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_nirspec_nrs1_wcs(_bigdata):
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    input_file = os.path.join(_bigdata, 'nirspec', 'test_wcs', 'nrs1-fs', 'jw00023001001_01101_00001_NRS1_ramp_fit.fits')
    ref_file = os.path.join(_bigdata, 'nirspec', 'test_wcs', 'nrs1-fs', 'jw00023001001_01101_00001_NRS1_ramp_fit_assign_wcs.fits')

    result = AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')
    result.close()

    im = ImageModel(result.meta.filename)
    imref = ImageModel(ref_file)

    for slit in ['S200A1', 'S200A2', 'S400A1', 'S1600A1']:
        w = nirspec.nrs_wcs_set_input(im, slit)
        grid = grid_from_bounding_box(w.bounding_box)
        ra, dec, lam = w(*grid)
        wref = nirspec.nrs_wcs_set_input(imref, slit)
        raref, decref, lamref = wref(*grid)

        assert_allclose(ra, raref, equal_nan=True)
        assert_allclose(dec, decref, equal_nan=True)
        assert_allclose(lam, lamref, equal_nan=True)
