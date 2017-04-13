from __future__ import absolute_import

import os
import pytest

from . import helpers
from .helpers import full_pool_rules

from .. import (generate, load_asn)

# Temporarily skip if running under Travis
pytestmark = pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason='Temporarily disable due to performance issues'
)


def test_generate(full_pool_rules):
    pool, rules, pool_fname = full_pool_rules
    (asns, orphaned) = generate(pool, rules)
    assert len(asns) == 11
    assert len(orphaned) == 2
    for asn in asns:
        asn_name, asn_store = asn.dump()
        asn_table = load_asn(asn_store)
        schemas = rules.validate(asn_table)
        assert len(schemas) > 0


def test_serialize(full_pool_rules):
    pool, rules, pool_fname = full_pool_rules
    (asns, orphaned) = generate(pool, rules)
    for asn in asns:
        for format in asn.ioregistry:
            fname, serialized = asn.dump(format=format)
            assert serialized is not None
            recovered = load_asn(serialized)
            assert recovered is not None


def test_unserialize():
    asn_file = helpers.t_path(
        'data/asn_mosaic.json'
    )
    with open(asn_file, 'r') as asn_fp:
        asn = load_asn(asn_fp)
    assert isinstance(asn, dict)
