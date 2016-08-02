from __future__ import absolute_import

from . import helpers
from .helpers import full_pool_rules

from .. import (Association, generate)
from ..association import SERIALIZATION_PROTOCOLS, validate


def test_generate(full_pool_rules):
    pool, rules = full_pool_rules
    (asns, orphaned) = generate(pool, rules)
    assert len(asns) == 4
    assert len(orphaned) == 0
    for asn in asns:
        asn_name, asn_store = asn.dump()
        asn_table = Association.load(asn_store)
        schemas = rules.validate(asn_table)
        assert len(schemas) > 0
        schemas = validate(asn_table)
        assert len(schemas) > 0


def test_serialize(full_pool_rules):
    pool, rules = full_pool_rules
    (asns, orphaned) = generate(pool, rules)
    for protocol in SERIALIZATION_PROTOCOLS:
        for asn in asns:
            fname, serialized = asn.dump(protocol=protocol)
            assert serialized is not None
            recovered = Association.load(serialized)
            assert recovered is not None


def test_unserialize():
    asn_file = helpers.t_path(
        'data/asn_mosaic.json'
    )
    with open(asn_file, 'r') as asn_fp:
        asn = Association.load(asn_fp)
    assert isinstance(asn, dict)
