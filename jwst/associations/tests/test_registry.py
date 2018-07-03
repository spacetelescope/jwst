"""Test the Registry"""
from .helpers import full_pool_rules

from .. import Association

def test_registry_match(full_pool_rules):
    """Test the match method"""
    pool, rules, pool_fname = full_pool_rules

    assert len(rules.schemas) > 0
    matches = rules.match(pool[1])
    assert isinstance(matches, tuple)
    asns = matches[0]
    reprocess_list = matches[1]
    assert isinstance(asns, list)
    assert isinstance(reprocess_list, list)
    assert len(asns) >= 1
