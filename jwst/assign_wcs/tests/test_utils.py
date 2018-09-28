"""
Test the utility functions
"""

import os

from astropy.table import QTable

from ...lib.catalog_utils import SkyObject

from ..util import get_object_info
from . import data

data_path = os.path.split(os.path.abspath(data.__file__))[0]


def get_file_path(filename):
    """
    Construct an absolute path.
    """
    return os.path.join(data_path, filename)


def read_catalog(catalogname): 
    return get_object_info(catalogname)


def test_create_grism_objects():
    source_catalog = get_file_path('step_SourceCatalogStep_cat.ecsv')

    # create from test ascii file
    grism_objects = read_catalog(source_catalog)
    assert isinstance(grism_objects, list), "return grism objects were not a list"

    required_fields = list(SkyObject()._fields)
    go_fields = grism_objects[0]._fields
    assert all([a == b for a, b in zip(required_fields, go_fields)]), "Required fields mismatch for SkyObject and GrismObject"

    # create from QTable object
    tempcat = QTable.read(source_catalog, format='ascii.ecsv')
    grism_object_from_table = read_catalog(tempcat)
    assert isinstance(grism_object_from_table, list), "return grism objects were not a list"
