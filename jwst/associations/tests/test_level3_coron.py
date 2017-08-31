"""Test Coronograpy rules

Coronograpy rules produce associations meant to be
processed by `calwebb_coron3`
"""
from __future__ import absolute_import

from .helpers import (
    BasePoolRule,
    PoolParams,
    combine_pools,
    registry_level3_only,
    t_path
)

from ..generate import generate


class TestLevel3Coron(BasePoolRule):

    pools = [
        PoolParams(
            path=t_path('data/pool_013_coron_nircam.csv'),
            n_asns=4,
            n_orphaned=2
        ),
    ]

    valid_rules = [
        'Asn_Coron',
    ]


def test_asn_type():
    """Ensure type is `coron3`"""

    pool = combine_pools(t_path('data/pool_013_coron_nircam.csv'))
    (asns, orphaned) = generate(pool, registry_level3_only())
    assert len(asns) > 0
    asn = asns[0]
    assert asn['asn_type'] == 'coron3'
