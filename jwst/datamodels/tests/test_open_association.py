import json
import os
import os.path
import warnings


from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer


# Define artificial memory size
MEMORY = 100  # 100 bytes
DATADIR = "data"

# Utilities
def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.join(os.path.dirname(__file__), DATADIR)
    return os.path.join(test_dir, partial_path)


def test_open_association():
    """Test for opening an association"""

    asn_file = t_path('association.json')
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(asn_file) as c:
            assert isinstance(c, ModelContainer)
            for model in c:
                assert model.meta.asn.table_name == "association.json"
                assert model.meta.asn.pool_name == "pool"


def test_container_open_asn_with_sourcecat():
    path = t_path("association_w_cat.json")
    with datamodels.open(path, asn_exptypes="science") as c:
        for model in c:
            assert model.meta.asn.table_name == "association_w_cat.json"


def test_open_with_relative_path_inside_asn():
    """Coverage for bug where relative paths inside filenames in an ASN would not be found,
    see JP-2038 / GitHub Issue 5950"""
    asn_file = t_path("association_with_paths.json")

    # ensure that there are indeed relative paths in the new asn
    with open(asn_file, "r") as f:
        asn_data = json.load(f)
        assert asn_data["products"][0]["members"][0]["expname"].split("/") == [DATADIR, "test.fits"]

    # cehck that this can be opened
    with datamodels.open(asn_file) as c:
        for model in c:
            assert model.meta.asn.table_name == "association_with_paths.json"
            assert model.meta.asn.pool_name == "pool"