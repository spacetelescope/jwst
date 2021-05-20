"""Test Level2 background nods"""
from jwst.associations.tests.helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from jwst.associations import generate
from jwst.associations.main import constrain_on_candidates

DITHER_PATTERN_MULTIPLIER = {
    '0': 1,  # No pattern, 1-to-1 exposure count
    '2': 2,  # Spatial, 2-to-1 exposure count
    '3': 3,  # Spectral, 3-to-1 exposure count
    '4': 4,  # Both, 4-to-1 exposure count
}

def test_nrs_msa_nod():
    pool = combine_pools(t_path('data/pool_023_nirspec_msa_3nod.csv'))
    all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(global_constraints=all_candidates))
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
    assert len(asns) == 30
    for asn in asns:
        n_dithers = int(asn.constraints['nods'].value)
        n_spectral_dithers = int(asn.constraints['subpxpts'].value)
        #  Expect self + all exposures not at the same primary dither
        n_members = n_dithers - n_spectral_dithers + 1
        assert len(asn['products'][0]['members']) == n_members


def test_nrs_fixedslit_nod_chop():
    """Test NIRSpec Fixed-slit background nods"""
    pool = combine_pools(t_path('data/pool_025_nirspec_fss_nod_chop.csv'))
    constraint_all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(
        global_constraints=constraint_all_candidates)
    )
    assert len(asns) == 8
    for asn in asns:
        assert asn['asn_rule'] in ['Asn_Lv2NRSFSS', 'Asn_Lv2SpecSpecial']
        if asn['asn_rule'] == 'Asn_Lv2SpecSpecial':
            assert len(asn['products'][0]['members']) == 1
        else:
            nods = int(asn.constraints['nods'].value)
            if asn['asn_id'].startswith('c'):
                nods += 1
            assert len(asn['products'][0]['members']) == nods
