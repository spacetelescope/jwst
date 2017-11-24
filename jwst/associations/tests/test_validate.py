import pytest

from . import helpers

from .. import (
    AssociationRegistry,
    AssociationNotValidError,
    load_asn
)


def test_invalid():
    rules = AssociationRegistry()
    with pytest.raises(AssociationNotValidError):
        rules.validate({})


def test_valid():
    rules = AssociationRegistry()
    asn_file = helpers.t_path(
        'data/test_image_asn.json'
    )
    with open(asn_file, 'r') as asn_fp:
        asn = load_asn(asn_fp)
    valid_schema_list = rules.validate(asn)
    assert isinstance(valid_schema_list, list)
