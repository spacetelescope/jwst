"""Test Level2 background nods"""
from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from .. import generate
from ..main import constrain_on_candidates


def test_nrs_msa_nod():
    pool = combine_pools(t_path('data/pool_023_nirspec_msa_3nod.csv'))
    all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(global_constraints=all_candidates))
    assert len(asns) == 12
    for asn in asns:
        assert len(asn['products'][0]['members']) == 3
