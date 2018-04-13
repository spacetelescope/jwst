import os
import pytest

import numpy as  np
from numpy.testing import utils
from asdf import AsdfFile
from astropy.io import fits

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel, fits_support

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_miri_ifu_wcs():
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    try:
        os.remove("miri_ifu_wcs_output.fits")
    except:
        pass

    input_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'ifu', 'jw00024001001_01101_00001_MIRIFUSHORT_uncal_MiriSloperPipeline.fits')
    ref_file = os.path.join(_bigdata, 'miri', 'test_wcs', 'ifu', 'jw00024001001_01101_00001_MIRIFUSHORT_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file='miri_ifu_wcs_output.fits'
                       )
    im = ImageModel('miri_ifu_wcs_output.fits')
    imref = ImageModel(ref_file)

    # Get the valid region
    f = fits.open('miri_ifu_wcs_output.fits')
    region = AsdfFile.open('/grp/crds/cache/references/jwst/{}'.format(f[0].header['R_REGION'].replace('crds://', '')))
    y, x = np.nonzero(region.tree['regions'])

    # Get indices where pixels == 0. These should be NaNs in the output.
    ind_zeros = region.tree['regions'] == 0

    ra, dec, lam = im.meta.wcs(x, y)
    raref, decref, lamref = imref.meta.wcs(x, y)
    utils.assert_allclose(ra, raref, equal_nan=True)
    utils.assert_allclose(dec, decref, equal_nan=True)
    utils.assert_allclose(lam, lamref, equal_nan=True)

    # Test that we got NaNs at ind_zero
    # DISABLED - What is this `all()` business?
    #assert(np.isnan(ra).nonzero()[0] == ind0[0])all()
    #assert(np.isnan(ra).nonzero()[1] == ind0[1])all()

    # Test the inverse transform
    x1, y1 = im.meta.wcs.backward_transform(ra, dec, lam)
    # DISABLED - What is this `all()` business?
    #assert(np.isnan(x1).nonzero()[0] == ind0[0])all()
    #assert (np.isnan(x1).nonzero()[1] == ind0[1])all()

    # Also run a smoke test with values outside the region.
    dec[100][200] = -80
    ra[100][200] = 7
    lam[100][200] = 15

    x2, y2 = im.meta.wcs.backward_transform(ra, dec, lam)
    assert np.isnan(x2[100][200])
    assert np.isnan(x2[100][200])

