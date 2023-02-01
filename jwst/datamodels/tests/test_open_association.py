import os
import os.path
import warnings


from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer


# Define artificial memory size
MEMORY = 100  # 100 bytes


# Utilities
def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.join(os.path.dirname(__file__), 'data')
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
