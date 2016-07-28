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
            path=t_path('data/pool_of_everything.csv'),
            n_asns=3,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_Image',
    ]
