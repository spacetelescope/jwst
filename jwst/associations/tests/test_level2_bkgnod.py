"""test_level3_dithers: Test of dither rules."""
from __future__ import absolute_import

from .helpers import BasePoolRule, PoolParams


class TestLevel2Bkgnod(BasePoolRule):

    pools = [
        PoolParams(
            path='tests/data/jw95070_20150615T143412_pool.csv',
            n_asns=2,
            n_orphaned=2,
            kwargs={'delimiter': ','}
        ),
    ]

    valid_rules = [
        'Asn_MIRI_LRS_BKGNOD',
    ]
