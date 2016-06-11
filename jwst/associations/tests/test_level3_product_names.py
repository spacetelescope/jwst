from __future__ import absolute_import

import re

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)

LEVEL3_PRODUCT_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '(?:-(?P<acid>[a-z]\d{3,4}))?'
    '_(?P<targname>(?:t\d{3})|(?:s\d{5}))'
    '(?:-(?P<epoch>epoch\d+))?'
    '_(?P<instrument>.+)'
    '_(?P<opt_elem>.+)'
    '_(?P<product_type>.+)'
    '\.fits'
)


class TestProductNames():
    pool_files = [
        helpers.t_path('data/jw93060_20150312T160130_pool.csv'),
        helpers.t_path('data/jw96090_20160406T233447_pool.csv')
    ]

    global_constraints = {
        'asn_candidate_ids': {
            'value': ['2'],
            'inputs': ['ASN_CANDIDATE_ID', 'OBS_NUM'],
            'force_unique': True,
        }
    }

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_level35_names(self):
        rules = AssociationRegistry()
        for pool_file in self.pool_files:
            pool = AssociationPool.read(pool_file)
            (asns, orphaned) = generate(pool, rules)
            for asn in asns:
                product_name = asn.data['products'][0]['name']
                m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
                yield helpers.not_none, m

    def test_level3_names(self):
        rules = AssociationRegistry(
            global_constraints=self.global_constraints
        )
        pool = AssociationPool.read(self.pool_files[0])
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            product_name = asn.data['products'][0]['name']
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            yield helpers.not_none, m
            yield helpers.check_equal, m.groupdict()['acid'], 'o002'

    def test_multiple_optelems(self):
        rules = AssociationRegistry()
        pool = AssociationPool.read(self.pool_files[1])
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            product_name = asn.data['products'][0]['name']
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            yield helpers.not_none, m
            if asn.constraints['opt_elem2']['value'] == 'CLEAR':
                yield helpers.check_not_in_list, '-', m.groupdict()['opt_elem']
            else:
                yield helpers.check_in_list, '-', m.groupdict()['opt_elem']
