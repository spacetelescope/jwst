"""Test general Level 3 rules environment"""
from __future__ import absolute_import

import re

from . import helpers
from .helpers import (combine_pools, t_path)

from .. import (AssociationRegistry, generate)


# Level 3 product name templates
L3_PRODUCT_NAME = (
    'jw(?P<program>\d{5})'
    '-(?P<acid>[a-z]\d{3,4})'
    '_(?P<targetid>t\d{3})'
    '_(?P<instrument>.+?)'
    '_(?P<opt_elem>.+?)'
    '_(?P<ptype>.+)\.fits'
)


class TestLevel3Environment(object):

    def test_l35_productname(self):
        rules = AssociationRegistry()
        pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
        asns, orphaned = generate(pool, rules)
        asn = asns[0]
        match = re.match(L3_PRODUCT_NAME, asn['products'][0]['name'])
        yield helpers.not_none, match
        matches = match.groupdict()
        yield helpers.check_equal, matches['program'], '99009'
        yield helpers.check_equal, matches['acid'], 'a3001'
        yield helpers.check_equal, matches['targetid'], 't001'
        yield helpers.check_equal, matches['instrument'], 'miri'
        yield helpers.check_equal, matches['opt_elem'], 'f560w'
        yield helpers.check_equal, matches['ptype'], '{product_type}'

    def test_l3_productname(self):
        global_constraints = {}
        global_constraints['asn_candidate_ids'] = {
            'value': '.+o001.+',
            'inputs': ['ASN_CANDIDATE'],
            'force_unique': True,
            'is_acid': True,
        }
        rules = AssociationRegistry(global_constraints=global_constraints)
        pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
        asns, orphaned = generate(pool, rules)
        asn = asns[0]
        match = re.match(L3_PRODUCT_NAME, asn['products'][0]['name'])
        yield helpers.not_none, match
        matches = match.groupdict()
        yield helpers.check_equal, matches['program'], '99009'
        yield helpers.check_equal, matches['acid'], 'o001'
        yield helpers.check_equal, matches['targetid'], 't001'
        yield helpers.check_equal, matches['instrument'], 'miri'
        yield helpers.check_equal, matches['opt_elem'], 'f560w'
        yield helpers.check_equal, matches['ptype'], '{product_type}'
