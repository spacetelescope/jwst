"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import
import pytest

from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)
from .. import generate


@pytest.mark.parametrize(
    'pool_name, n_asns',
    [
        (t_path('data/pool_011_spec_miri_lv2bkg_lrs.csv'), 3),
        (t_path('data/pool_009_spec_miri_lv2bkg.csv'), 16),
        (t_path('data/pool_010_spec_nirspec_lv2bkg.csv'), 12),
    ]
)
def test_level2_bkg_cand(pool_name, n_asns):
    rules = registry_level2_only()
    pool = combine_pools(pool_name)
    asns, orphaned = generate(pool, rules)
    assert len(asns) == n_asns
    for asn in asns:
        assert asn['asn_type'] == 'spec2'
        assert asn['asn_rule'] in ['Asn_Lv2Spec', 'Asn_Lv2SpecBkg']

        members = asn['members']
        if len(members) == 1:
            assert members[0]['exptype'] == 'SCIENCE'
        else:
            exptypes = [
                member['exptype']
                for member in members
            ]
            assert exptypes.count('SCIENCE') == 1
            assert exptypes.count('BACKGROUND') >= 1


@pytest.mark.xfail(reason='No determined yet')
def test_level2_bkgnod():
    assert False

"""
class TestLevel2Bkgnod(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/jw95070_20150615T143412_pool.csv'),
            n_asns=2,
            n_orphaned=2,
            kwargs={'delimiter': ','}
        ),
    ]

    valid_rules = [
        'Asn_MIRI_LRS_BKGNOD',
    ]
"""
