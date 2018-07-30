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
    constraint_all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(
        global_constraints=constraint_all_candidates)
    )
    assert len(asns) == 12
    for asn in asns:
        assert len(asn['products'][0]['members']) == 3


def test_nrs_fixedslit_nod():
    """Test NIRSpec Fixed-slit background nods"""
    pool = combine_pools(t_path('data/pool_024_nirspec_fss_nods.csv'))
    constraint_all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(
        global_constraints=constraint_all_candidates)
    )
    assert len(asns) == 6
    for asn in asns:
        nods = int(asn.constraints['nods'].value)
        assert len(asn['products'][0]['members']) == nods


def test_nrs_fixedslit_nod_chop():
    """Test NIRSpec Fixed-slit background nods"""
    pool = combine_pools(t_path('data/pool_025_nirspec_fss_nod_chop.csv'))
    constraint_all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(
        global_constraints=constraint_all_candidates)
    )
    assert len(asns) == 14
    for asn in asns:
        assert asn['asn_rule'] in ['Asn_Lv2NRSFSS', 'Asn_Lv2SpecSpecial']
        if asn['asn_rule'] == 'Asn_Lv2SpecSpecial':
            assert len(asn['products'][0]['members']) == 1
        else:
            nods = int(asn.constraints['nods'].value)
            if asn['asn_id'].startswith('c'):
                nods += 1
            assert len(asn['products'][0]['members']) == nods
