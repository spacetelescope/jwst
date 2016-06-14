"""Test general Level 3 rules environment"""
from __future__ import absolute_import

import re

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)


# Level 3 product name templates
L35_PRODUCT_NAME = (
    'jw(?P<program>\d{5})'
    '_(?P<targetid>t\d{3})'
    '_(?P<instrument>.+)'
    '_(?P<opt_elem>.+)'
    '_(?P<ptype>.+)\.fits'
)
L3_PRODUCT_NAME = (
    'jw(?P<program>\d{5})'
    '-(?P<acid>[a-z]\d{3,4})'
    '_(?P<targetid>t\d{3})'
    '_(?P<instrument>.+)'
    '_(?P<opt_elem>.+)'
    '_(?P<ptype>.+)\.fits'
)


class TestLevel3Environment(object):

    pools_size = [
        (
            helpers.t_path('data/jw93060_20150312T160130_pool.csv'),
            14
        ),
    ]

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_l35_productname(self):
        rules = AssociationRegistry()
        pool = AssociationPool.read(self.pools_size[0][0])
        asns, orphaned = generate(pool, rules)
        asn = asns[0]
        match = re.match(L35_PRODUCT_NAME, asn.data['products'][0]['name'])
        yield helpers.not_none, match
        matches = match.groupdict()
        yield helpers.check_equal, matches['program'], '93060'
        yield helpers.check_equal, matches['targetid'], 't001'
        yield helpers.check_equal, matches['instrument'], 'miri'
        yield helpers.check_equal, matches['opt_elem'], 'f560w'
        yield helpers.check_equal, matches['ptype'], 'dither'

    def test_l3_productname(self):
        global_constraints = {}
        global_constraints['asn_candidate_ids'] = {
            'value': '1',
            'inputs': ['ASN_CANDIDATE_ID', 'OBS_NUM'],
            'force_unique': True,
        }
        rules = AssociationRegistry(global_constraints=global_constraints)
        pool = AssociationPool.read(self.pools_size[0][0])
        asns, orphaned = generate(pool, rules)
        asn = asns[0]
        match = re.match(L3_PRODUCT_NAME, asn.data['products'][0]['name'])
        yield helpers.not_none, match
        matches = match.groupdict()
        yield helpers.check_equal, matches['program'], '93060'
        yield helpers.check_equal, matches['acid'], 'o001'
        yield helpers.check_equal, matches['targetid'], 't001'
        yield helpers.check_equal, matches['instrument'], 'miri'
        yield helpers.check_equal, matches['opt_elem'], 'f560w'
        yield helpers.check_equal, matches['ptype'], 'dither'
