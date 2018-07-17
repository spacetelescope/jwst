from . import helpers
from .helpers import full_pool_rules

from .. import (
    AssociationPool,
    AssociationRegistry,
    generate,
    load_asn
)


def test_simple():
    """Test generate on simple registry"""
    registry = AssociationRegistry(
        [helpers.t_path('data/rules_basic.py')],
        include_default=False
    )
    pool = AssociationPool()
    pool['value'] = ['row1', 'row2']

    asns = generate(pool, registry)
    assert len(asns) == 1
    assert len(asns[0]['members']) == 2


@helpers.runslow
def test_generate(full_pool_rules):
    pool, rules, pool_fname = full_pool_rules
    asns = generate(pool, rules)
    assert len(asns) == 37
    for asn in asns:
        asn_name, asn_store = asn.dump()
        asn_table = load_asn(asn_store)
        schemas = rules.validate(asn_table)
        assert len(schemas) > 0


@helpers.runslow
def test_serialize(full_pool_rules):
    pool, rules, pool_fname = full_pool_rules
    asns = generate(pool, rules)
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
