import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box
from ci_watson.artifactory_helpers import get_bigdata

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel


testdata = [
    ('nrs1', 'jw00011001001_01120_00001_NRS1_rate.fits', 'jw00011001001_01120_00001_NRS1_assign_wcs.fits'),
    ('nrs1_opaque', 'jw00011001001_01120_00001_NRS1_rate_opaque.fits', 'jw00011001001_01120_00001_NRS1_rate_opaque_assign_wcs.fits'),
    ('nrs2', 'NRSIFU-COMBO-030_NRS2_SloperPipeline.fits', 'NRSIFU-COMBO-030_NRS2_SloperPipeline_assign_wcs.fits')
]

@pytest.mark.bigdata
@pytest.mark.parametrize("test_id, input_file, truth_file", testdata)
def test_nirspec_ifu_wcs(envopt, _jail, test_id, input_file, truth_file):
    """
    Regression test of creating a WCS object and doing pixel to sky transformation.
    """
    del test_id

    input_file = get_bigdata('jwst-pipeline', envopt,
                             'nirspec', 'test_wcs', 'nrs1-ifu', input_file)
    truth_file = get_bigdata('jwst-pipeline', envopt,
                             'nirspec', 'test_wcs', 'nrs1-ifu', 'truth', truth_file)

    result = AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')
    result.close()

    im = ImageModel(result.meta.filename)
    imref = ImageModel(truth_file)
    w = nirspec.nrs_wcs_set_input(im, 0)
    grid = grid_from_bounding_box(w.bounding_box)
    ra, dec, lam = w(*grid)
    wref = nirspec.nrs_wcs_set_input(imref, 0)
    raref, decref, lamref = wref(*grid)

    # equal_nan is used here as many of the entries are nan.
    # The domain is defined but it is only a few entries in there that are valid
    # as it is a curved narrow slit.
    assert_allclose(ra, raref, equal_nan=True)
    assert_allclose(dec, decref, equal_nan=True)
    assert_allclose(lam, lamref, equal_nan=True)
