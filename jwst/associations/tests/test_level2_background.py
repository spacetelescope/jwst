"""Test Level2 background nods"""
from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from .. import generate


def test_nrs_msa_nod():
    pool = combine_pools(t_path('data/pool_023_nirspec_msa_3nod.csv'))
    asns = generate(pool, registry_level2_only())
    assert len(asns) == 12
    for asn in asns:
        assert len(asn['products'][0]['members']) == 3
