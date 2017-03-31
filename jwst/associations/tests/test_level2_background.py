"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import
import pytest

from .helpers import (
    combine_pools,
    registry_level2_only,
    t_path
)
from .. import generate


def test_level2_background_candidate():
    rules = registry_level2_only()
    pool = combine_pools(
        t_path(
            'data/pool_009_spec_miri_lrs_fsbkg.csv'
        )
    )
    asns, orphaned = generate(pool, rules)
    assert len(asns) == 1
    members = asns[0]['members']
    assert len(members) == 2
    n_bkgs = 0
    for member in members:
        if 'bkgexps' in members:
            n_bkgs += 1
    assert n_bkgs == 1


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
