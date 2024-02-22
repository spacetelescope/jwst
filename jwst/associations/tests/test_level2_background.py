"""Test Level2 background nods"""
from jwst.associations.tests.helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)

from jwst.associations import generate
from jwst.associations.lib.utilities import constrain_on_candidates

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


def test_nrs_msa_nod_subpix():
    # This is a pool with NIRSpec MOS exposures at 3 primary nods and 3 subpixel
    # dither positions at each nod, for a total 9 exposures on each of two
    # (nrs1, nrs2) detectors, so a grand total of 18 exposures and 18 asns
    # (one spec2 asn for each exposure). So each spec2 asn should contain 1
    # science member, and the 6 exposures from the other 2 primary nod positions
    # as background members.
    pool = combine_pools(t_path('data/pool_023b_nirspec_msa_3nod_subpix.csv'))
    all_candidates = constrain_on_candidates(None)
    asns = generate(pool, registry_level2_only(global_constraints=all_candidates))
    assert len(asns) == 18
    for asn in asns:
        assert len(asn['products'][0]['members']) == 7
        nsci = 0
        nbkg = 0
        for member in asn['products'][0]['members']:
            if member['exptype'] == 'science':
                nsci += 1
            if member['exptype'] == 'background':
                nbkg += 1
        assert nsci == 1
        assert nbkg == 6


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
