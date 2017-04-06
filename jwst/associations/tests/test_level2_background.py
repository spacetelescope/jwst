"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import
import pytest

from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)
from .. import generate


def test_level2_bkg_cand_miri():
    rules = registry_level2_only()
    pool = combine_pools(
        t_path(
            'data/pool_009_spec_miri_lv2bkg.csv'
        )
    )
    asns, orphaned = generate(pool, rules)
    assert len(asns) == 2
    for asn in asns:
        members = asn['members']
        if asn['asn_type'] == 'Asn_Lv2Spec':
            assert len(members) == 1
        elif asn['asn_type'] == 'AsnLv2SpecBkg':
            assert len(members) == 2
            exptypes = set(
                member['exptype']
                for member in members
            )
            assert exptypes.issuperset('SCIENCE', 'BACKGROUND')


def test_level2_bkg_cand_nrs():
    rules = registry_level2_only()
    pool = combine_pools(
        t_path(
            'data/pool_010_spec_nirspec_lv2bkg.csv'
        )
    )
    asns, orphaned = generate(pool, rules)
    assert len(asns) == 6
    for asn in asns:
        members = asn['members']
        if asn['asn_type'] == 'Asn_Lv2Spec':
            assert len(members) >= 2
        elif asn['asn_type'] == 'AsnLv2SpecBkg':
            assert len(members) == 2
            exptypes = set(
                member['exptype']
                for member in members
            )
            assert exptypes.issuperset('SCIENCE', 'BACKGROUND')


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
