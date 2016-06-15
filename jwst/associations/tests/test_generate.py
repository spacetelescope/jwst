from __future__ import absolute_import

import pytest

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)
from ..association import SERIALIZATION_PROTOCOLS


@pytest.fixture
def pool_rules():
    pool_file = helpers.t_path('data/jw93060_20150312T160130_pool.csv')
    rules = AssociationRegistry()
    pool = AssociationPool.read(pool_file)
    return (pool, rules)


def test_generate(pool_rules):
    pool, rules = pool_rules
    (asns, orphaned) = generate(pool, rules)
    assert len(asns) == 14
    assert len(orphaned) == 36


def test_serialize():
    pool, rules = pool_rules()
    (asns, orphaned) = generate(pool, rules)
    for protocol in SERIALIZATION_PROTOCOLS:
        for asn in asns:
            yield helpers.not_none, asn.serialize(protocol=protocol)
