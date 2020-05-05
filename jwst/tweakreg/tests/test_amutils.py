"""Test astrometric utility functions for alignment"""
import pytest

import asdf
import numpy as np

from jwst.tweakreg import astrometric_utils as amutils

# Define input GWCS specification to be used for these tests
WCS_NAME = 'mosaic_long_i2d_gwcs.asdf'  # Derived using B7.5 Level 3 product
EXPECTED_NUM_SOURCES = 2347
EXPECTED_RADIUS = 0.02564497890604383
TEST_CATALOG = 'GAIADR2'

# Set test sensitivity
TEST_RTOL = 1e-6

# Utility functions for the tests
def read_test_gwcs():
    asdf_file = asdf.open(WCS_NAME)
    wcsobj = asdf_file['wcs']

    return wcsobj


# Test definitions
def test_radius():
    # Read GWCS object
    wcsobj = read_test_gwcs()
    
    # compute radius
    radius, fiducial = amutils.compute_radius(wcsobj)
    
    # check results
    assert(np.allclose([radius], [EXPECTED_RADIUS], rtol=TEST_RTOL))
    
    
def test_get_catalog():
    # Read GWCS object
    wcsobj = read_test_gwcs()

    # Interpret the GWCS 
    radius, fiducial = amutils.compute_radius(wcsobj)
    
    # Get the catalog
    cat = amutils.get_catalog(fiducial[0], fiducial[1], 
                              sr=radius, 
                              catalog=TEST_CATALOG)
    # check results
    assert(len(cat) == EXPECTED_NUM_SOURCES)


def test_create_catalog():
    # Read GWCS object 
    wcsobj = read_test_gwcs()
    
    # Create catalog
    gcat = amutils.create_astrometric_catalog(None, 
                                                existing_wcs=wcsobj,
                                                catalog=TEST_CATALOG,
                                                output=None)    

    # check that we got expected number of sources
    assert(len(gcat) == EXPECTED_NUM_SOURCES)
  
