"""Test Coronograpy rules

Coronograpy rules produce associations meant to be
processed by `calwebb_coron3`
"""
from __future__ import absolute_import

from .helpers import BasePoolRule, PoolParams, t_path


class TestLevel3Coron(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/pool_013_coron_nircam.csv'),
            n_asns=2,
            n_orphaned=1
        ),
    ]

    valid_rules = [
        'Asn_Coron',
    ]
