"""Test general asn_generate operations"""
import pytest

from jwst.associations import (
    generate,
    load_asn
)


@pytest.mark.slow
def test_generate(full_pool_rules):
    """Run a full sized pool using all rules"""
    pool, rules, pool_fname = full_pool_rules
    asns = generate(pool, rules)
    assert len(asns) == 97
    for asn in asns:
        asn_name, asn_store = asn.dump()
        asn_table = load_asn(asn_store)
        schemas = rules.validate(asn_table)
        assert len(schemas) > 0


@pytest.mark.slow
def test_serialize(full_pool_rules):
    """Test serializing roundtripping"""
    pool, rules, pool_fname = full_pool_rules
    asns = generate(pool, rules)
    for asn in asns:
        for format in asn.ioregistry:
            fname, serialized = asn.dump(format=format)
            assert serialized is not None
            recovered = load_asn(serialized)
            assert recovered is not None
