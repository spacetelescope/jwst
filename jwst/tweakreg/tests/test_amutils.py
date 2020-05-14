"""Test astrometric utility functions for alignment"""
import os

import asdf
import numpy as np
import pytest

from jwst.tweakreg import astrometric_utils as amutils


# Define input GWCS specification to be used for these tests
WCS_NAME = 'mosaic_long_i2d_gwcs.asdf'  # Derived using B7.5 Level 3 product
EXPECTED_NUM_SOURCES = 2347
EXPECTED_RADIUS = 0.02564497890604383
TEST_CATALOG = 'GAIADR2'


@pytest.fixture(scope="module")
def wcsobj():
    path = os.path.join(os.path.dirname(__file__), WCS_NAME)
    with asdf.open(path) as asdf_file:
        wcs = asdf_file['wcs']

    return wcs


def test_radius(wcsobj):
    # compute radius
    radius, fiducial = amutils.compute_radius(wcsobj)

    # check results
    np.testing.assert_allclose(radius, EXPECTED_RADIUS, rtol=1e-6)


def test_get_catalog(wcsobj):
    # Get radius and fiducial
    radius, fiducial = amutils.compute_radius(wcsobj)

    # Get the catalog
    cat = amutils.get_catalog(fiducial[0], fiducial[1], sr=radius,
                              catalog=TEST_CATALOG)

    assert len(cat) == EXPECTED_NUM_SOURCES


def test_create_catalog(wcsobj):
    # Create catalog
    gcat = amutils.create_astrometric_catalog(None, existing_wcs=wcsobj,
        catalog=TEST_CATALOG, output=None)

    # check that we got expected number of sources
    assert len(gcat) == EXPECTED_NUM_SOURCES
