"""Test general Level 3 rules environment"""
import pytest

from .helpers import (
    combine_pools,
    registry_level3_only,
    t_path
)

from .. import generate


def test_meta():
    rules = registry_level3_only()
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns = generate(pool, rules)
    assert len(asns) == 1
    asn = asns[0]
    data = asn.data
    assert data['program'] == '99009'
    assert data['target'] == 't001'
    assert data['asn_type'] == 'image3'
    assert data['asn_id'] == 'a3001'
    assert data['asn_pool'] == 'pool_002_image_miri'
    assert data['asn_rule'] == 'Asn_Lv3Image'
    assert data['degraded_status'] == 'No known degraded exposures in association.'
    assert data['version_id'] is None
    assert data['constraints'] is not None


@pytest.mark.parametrize(
    'pool_file',
    [
        'data/pool_005_spec_niriss.csv',
        'data/pool_006_spec_nirspec.csv',
        'data/pool_007_spec_miri.csv',
        'data/pool_010_spec_nirspec_lv2bkg.csv',
        'data/pool_015_spec_nirspec_lv2bkg_reversed.csv',
        'data/pool_016_spec_nirspec_lv2bkg_double.csv',
        'data/pool_017_spec_nirspec_lv2imprint.csv',
    ]
)
def test_targacq(pool_file):
    """Test for existence of target acquisitions in associatons"""
    rules = registry_level3_only()
    pool = combine_pools(t_path(pool_file))
    asns = generate(pool, rules)
    assert len(asns) > 0
    for asn in asns:
        # Ignore reprocessed asn's with only science
        if not asn['asn_rule'] in ["Asn_Lv3SpecAux", "Asn_Lv3NRSIFUBackground"]:
            for product in asn['products']:
                exptypes = [
                    member['exptype'].lower()
                    for member in product['members']
                    ]
                assert 'target_acquisition' in exptypes
