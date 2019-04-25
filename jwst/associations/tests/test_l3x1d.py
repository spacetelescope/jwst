"""test_level3 MIR IFU rules"""
import pytest
import re

from .helpers import (
    BasePoolRule,
    PoolParams,
    func_fixture,
    generate_params,
    combine_pools,
    registry_level3_only,
    t_path
)

from .. import (AssociationPool, generate)
from ..lib.dms_base import DMSAttrConstraint


class TestL3x1d():
    def setUp(self):
        pass

    def tearDown(self):
        pass

#There where 9 associations and 10 orphaned items found.

pool_file = func_fixture(
    generate_params,
    scope='module',
    params=[
        t_path('data/jw00623_pool.csv'),
    ]
)


global_constraints = func_fixture(
    generate_params,
    scope='module',
    params=[
        DMSAttrConstraint(
            name='asn_candidate',
            value=['.+o002.+'],
            sources=['asn_candidate'],
            force_unique=True,
            is_acid=True,
            evaluate=True,
        ),
    ]
)

def test_l3_MIRFU():
    rules = registry_level3_only()
    pool = combine_pools(t_path('data/jw00623_pool.csv'))
    asns = generate(pool, rules)
