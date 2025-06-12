import pytest
from astropy.utils.data import get_pkg_data_filename

from jwst.associations import AssociationRegistry, AssociationNotValidError, load_asn


def test_invalid():
    rules = AssociationRegistry()
    with pytest.raises(AssociationNotValidError):
        rules.validate({})


def test_valid():
    rules = AssociationRegistry()
    asn_file = get_pkg_data_filename("data/test_image_asn.json", package="jwst.associations.tests")
    with open(asn_file, "r") as asn_fp:
        asn = load_asn(asn_fp)
    valid_schema_list = rules.validate(asn)
    assert isinstance(valid_schema_list, list)
