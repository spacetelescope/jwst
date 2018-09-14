import os
import pytest

import numpy as  np
from numpy.testing import utils

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel, RegionsModel
from jwst.stpipe import crds_client

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_miri_ifu_wcs(_bigdata):
    """
    Regression test of creating a WCS object and doing pixel to sky transformation.
    """
    input_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'ifu', 'jw00024001001_01101_00001_MIRIFUSHORT_uncal_MiriSloperPipeline.fits')
    ref_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'ifu', 'jw00024001001_01101_00001_MIRIFUSHORT_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file='miri_ifu_wcs', suffix='output'
                       )
    im = ImageModel('miri_ifu_wcs_output.fits')
    imref = ImageModel(ref_file)

    # Get the region file
    region = RegionsModel(crds_client.get_reference_file(im, 'regions'))
    
    # inputs
    shape = region.regions.shape
    y, x = np.mgrid[ : shape[0], : shape[1]]
    
    # Get indices where pixels == 0. These should be NaNs in the output.
    ind_zeros = region.regions == 0

    ra, dec, lam = im.meta.wcs(x, y)
    raref, decref, lamref = imref.meta.wcs(x, y)
    utils.assert_allclose(ra, raref, equal_nan=True)
    utils.assert_allclose(dec, decref, equal_nan=True)
    utils.assert_allclose(lam, lamref, equal_nan=True)

    # Test that we got NaNs at ind_zero
    assert(np.isnan(ra).nonzero()[0] == ind_zeros.nonzero()[0]).all()
    assert(np.isnan(ra).nonzero()[1] == ind_zeros.nonzero()[1]).all()

    # Test the inverse transform
    x1, y1 = im.meta.wcs.backward_transform(ra, dec, lam)
    assert(np.isnan(x1).nonzero()[0] == ind_zeros.nonzero()[0]).all()
    assert (np.isnan(x1).nonzero()[1] == ind_zeros.nonzero()[1]).all()

    # Also run a smoke test with values outside the region.
    dec[100][200] = -80
    ra[100][200] = 7
    lam[100][200] = 15

    x2, y2 = im.meta.wcs.backward_transform(ra, dec, lam)
    assert np.isnan(x2[100][200])
    assert np.isnan(x2[100][200])
