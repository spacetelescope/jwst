# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

import os
import shutil
import tempfile

import numpy as np
from numpy.testing import assert_array_almost_equal

from .. import ImageModel


TMP_DIR = None
FITS_FILE = None
TMP_FITS = None


def setup():
    global FITS_FILE, TMP_DIR, TMP_FITS, TMP_YAML, TMP_JSON
    ROOT_DIR = os.path.dirname(__file__)
    FITS_FILE = os.path.join(ROOT_DIR, 'sip.fits')

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')


def teardown():
    shutil.rmtree(TMP_DIR)


def _header_to_dict(x):
    return dict((a, b) for (a, b, c) in x)


def test_wcs():
    with ImageModel(FITS_FILE) as dm:
        wcs1 = dm.get_fits_wcs()
        dm2 = dm.copy()
        wcs2 = dm2.get_fits_wcs()

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2world(x, 1)
    world2 = wcs2.all_pix2world(x, 1)

    assert_array_almost_equal(world1, world2)

    wcs1.wcs.crpix[0] = 42.0

    dm2.set_fits_wcs(wcs1)
    header = _header_to_dict(dm2.extra_fits.PRIMARY.header)
    assert header['CRPIX1'] == 42.0

    wcs2 = dm2.get_fits_wcs()
    assert wcs2.wcs.crpix[0] == 42.0

    dm2.to_fits(TMP_FITS, clobber=True)

    with ImageModel(TMP_FITS) as dm3:
        wcs3 = dm3.get_fits_wcs()

    assert wcs3.wcs.crpix[0] == 42.0

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2world(x, 1)
    world2 = wcs3.all_pix2world(x, 1)

    dm4 = ImageModel()
    dm4.set_fits_wcs(wcs3)
    dm4.to_fits(TMP_FITS, clobber=True)

    with ImageModel(TMP_FITS) as dm5:
        wcs5 = dm5.get_fits_wcs()

    assert wcs5.wcs.crpix[0] == 42.0
