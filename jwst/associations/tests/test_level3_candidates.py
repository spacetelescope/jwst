"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import

from nose import SkipTest

from .helpers import BasePoolRule, PoolParams, t_path


class TestLevel3Candidates(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/jw96090_20160406T233447_pool.csv'),
            n_asns=2,
            n_orphaned=0,
            n_candidates=[1, 1],
        ),
    ]

    valid_rules = [
        'Asn_Image',
    ]
