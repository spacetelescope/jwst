import warnings

from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

# Define artificial memory size
MEMORY = 100  # 100 bytes


def test_open_association():
    """Test for opening an association"""

    asn_file = get_pkg_data_filename("data/association.json", package="jwst.datamodels.tests")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(asn_file) as c:
            assert isinstance(c, ModelContainer)
            for model in c:
                assert model.meta.asn.table_name == "association.json"
                assert model.meta.asn.pool_name == "pool"


def test_container_open_asn_with_sourcecat():
    path = get_pkg_data_filename("data/association_w_cat.json", package="jwst.datamodels.tests")
    with datamodels.open(path, asn_exptypes="science") as c:
        for model in c:
            assert model.meta.asn.table_name == "association_w_cat.json"
