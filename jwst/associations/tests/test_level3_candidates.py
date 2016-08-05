"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import

from .helpers import BasePoolRule, PoolParams, t_path


class TestLevel3Candidates(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/pool_001_candidates.csv'),
            n_asns=2,
            n_orphaned=0,
            candidates=[['100'], ['100']],
        ),
    ]

    valid_rules = [
        'Asn_Image',
    ]
