import os
import pytest
import numpy as np
from numpy.testing import utils

from jwst.assign_wcs import AssignWcsStep, nirspec

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]
from jwst.datamodels import ImageModel

from ..helpers import add_suffix


def test_nirspec_nrs1_wcs(_bigdata):
    """

    Regression test of creating a WCS object and doing pixel to sky transformation.

    """
    output_file_base, output_file = add_suffix('nrs1_fs_wcs_output.fits', 'assignwcsstep')

    try:
        os.remove(output_file)
    except:
        pass

    input_file = os.path.join(_bigdata, 'nirspec', 'test_wcs', 'nrs1-fs', 'jw00023001001_01101_00001_NRS1_ramp_fit.fits')
    ref_file = os.path.join(_bigdata, 'nirspec', 'test_wcs', 'nrs1-fs', 'jw00023001001_01101_00001_NRS1_assign_wcs.fits')

    AssignWcsStep.call(input_file,
                       output_file=output_file_base, name='assignwcsstep'
                       )
    im = ImageModel(output_file)
    imref = ImageModel(ref_file)
    #ystart = im.meta.subarray.ystart
    #yend = im.meta.subarray.ystart + im.meta.subarray.ysize-1
    #xstart = im.meta.subarray.xstart
    #xend = im.meta.subarray.xstart + im.meta.subarray.xsize -1
    #x, y = np.mgrid[ystart:yend, xstart: xend]
    for slit in ['S200A1', 'S200A2', 'S400A1', 'S1600A1']:
        w = nirspec.nrs_wcs_set_input(im, slit)
        #ra, dec, lam = getattr(im.meta, 'wcs_'+slit)(y, x)
        #raref, decref, lamref = getattr(imref.meta, 'wcs_'+slit)(y, x)
        y, x = np.mgrid[w.bounding_box[1][0]:w.bounding_box[1][1],
                        w.bounding_box[0][0]: w.bounding_box[0][1]]
        ra, dec, lam = w(x, y)
        wref = nirspec.nrs_wcs_set_input(im, slit)
        raref, decref, lamref = wref(x, y)
        utils.assert_allclose(ra, raref, equal_nan=True)
        utils.assert_allclose(dec, decref, equal_nan=True)
        utils.assert_allclose(lam, lamref, equal_nan=True)
