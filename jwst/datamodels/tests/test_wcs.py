import os
import warnings

import numpy as np
from numpy.testing import assert_array_almost_equal

from ..wcs_ref_models import FilteroffsetModel
from .. import ImageModel

import pytest


FITS_FILE = os.path.join(os.path.dirname(__file__), 'data', 'sip.fits')


def test_wcs(tmpdir):
    with ImageModel(FITS_FILE) as dm:

        # Refer to the data array to initialize it.
        dm.data = np.zeros((5, 5))

        # Now continue with the test.
        wcs1 = dm.get_fits_wcs()
        dm2 = dm.copy()
        wcs2 = dm2.get_fits_wcs()

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2world(x, 1)
    world2 = wcs2.all_pix2world(x, 1)

    assert_array_almost_equal(world1, world2)

    wcs1.wcs.crpix[0] = 42.0

    dm2.set_fits_wcs(wcs1)
    assert dm2.meta.wcsinfo.crpix1 == 42.0

    wcs2 = dm2.get_fits_wcs()
    assert wcs2.wcs.crpix[0] == 42.0

    dm2_tmp_fits = str(tmpdir.join("tmp_dm2.fits"))
    dm2.to_fits(dm2_tmp_fits)

    with ImageModel(dm2_tmp_fits) as dm3:
        wcs3 = dm3.get_fits_wcs()

    assert wcs3.wcs.crpix[0] == 42.0

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2world(x, 1)
    world2 = wcs3.all_pix2world(x, 1)

    dm4 = ImageModel()
    dm4.set_fits_wcs(wcs3)
    dm4_tmp_fits = str(tmpdir.join("tmp_dm4.fits"))
    dm4.to_fits(dm4_tmp_fits, overwrite=True)

    warnings.filterwarnings(action="ignore", message="The WCS transformation has more axes")
    with ImageModel(dm4_tmp_fits) as dm5:
        wcs5 = dm5.get_fits_wcs()

    assert wcs5.wcs.crpix[0] == 42.0


def test_wcs_ref_models():
    filters = [{'name': 'F090W', 'row_offset': 1, 'column_offset': 1},
               {'name': 'F070W', 'row_offset': 2, 'column_offset': 2}
              ]
    with FilteroffsetModel(filters=filters, instrument='NIRCAM', strict_validation=True) as fo:
        fo.filters == filters
        with pytest.raises(ValueError) as e:
            fo.validate()
            assert str(e.value()) == ("Model.meta is missing values for "
                                     "['description', 'reftype', 'author', 'pedigree',"
                                     "'useafter']")
                                     
    filters = [{'filter': 'F090W', 'pupil': 'GRISMR',
                'row_offset': 1, 'column_offset': 1},
               {'filter': 'F070W', 'pupil': 'GRISMC',
                'row_offset': 2, 'column_offset': 2}
              ]
    with FilteroffsetModel(filters=filters, instrument='NIRCAM', strict_validation=True) as fo:
        fo.filters == filters
        fo.meta.description = "Filter offsets"
        fo.meta.reftype = "filteroffset"
        fo.meta.author = "Unknown"
        fo.meta.pedigree = "GROUND"
        fo.meta.useafter = "2019-12-01"

        with pytest.raises(ValueError) as e:
            fo.validate()
            assert str(e.value()) == ('Expected detector for instrument NIRCAM to be one of "NRCA1, NRCA2, "'
                                      'NRCA3, NRCA4, NRCALONG, NRCB1, NRCB2, NRCB3, NRCB4, NRCBLONG"')
        fo.meta.instrument.detector = 'NRCA3'
        fo.validate()
