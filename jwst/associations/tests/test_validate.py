from __future__ import absolute_import

import pytest

from . import helpers

from ..association import (
    Association,
    AssociationRegistry,
    AssociationNotValidError
)


def test_invalid():
    rules = AssociationRegistry()
    with pytest.raises(AssociationNotValidError):
        rules.validate({})


def test_valid():
    rules = AssociationRegistry()
    asn_file = helpers.t_path(
        'data/jw96090_20160615t210324_mosaic_001_asn.json'
    )
    with open(asn_file, 'r') as asn_fp:
        asn = Association.load(asn_fp)
    valid_schema_list = rules.validate(asn)
    assert isinstance(valid_schema_list, list)
