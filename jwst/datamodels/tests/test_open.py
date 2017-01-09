"""
Test datamodel.open
"""

import os

from .. import (DataModel, ModelContainer)
from ..util import open


def test_open_fits():
    """Test opening a model from a FITS file"""

    fits_file = t_path('test.fits')
    m = open(fits_file)
    assert isinstance(m, DataModel)


def test_open_association():
    """Test for opening an association"""

    asn_file = t_path('association.json')
    m = open(asn_file)
    assert isinstance(m, ModelContainer)


# Utilities
def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)
