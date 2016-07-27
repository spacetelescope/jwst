"""Test Image rules

Image rules produce associations meant to be
processed by CALIMAGE3.

Such associations are mosaics and dithers.
"""
from __future__ import absolute_import

from .helpers import BasePoolRule, PoolParams, t_path


class TestLevel3Image(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/jw93065_20151105T004034_pool.csv'),
            n_asns=3,
            n_orphaned=0
        ),
        PoolParams(
            path=t_path('data/jw82600_001_20151107T165901_pool.csv'),
            n_asns=11,
            n_orphaned=298
        ),

    ]

    valid_rules = [
        'Asn_Image',
    ]
