from __future__ import absolute_import

import re

from . import helpers

from jwst.associations.association import AssociationRegistry
from jwst.associations.pool import AssociationPool
from jwst.associations.generate import generate

LEVEL3_ASN_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '(?:-(?P<acid>\d{4}))?'
    '_(?P<stamp>\d{8}(T|t)\d{6})'
    '_(?P<asn_type>.+)'
    '_(?P<sequence>\d{3})'
    '_asn'
)


class TestASNtNames():
    pool_file = 'tests/data/jw93060_20150312T160130_pool.csv'

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
        pool = AssociationPool.read(self.pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            asn_name = asn.asn_name
            m = re.match(LEVEL3_ASN_NAME_REGEX, asn_name)
            yield helpers.not_none, m

    def test_level3_names(self):
        rules = AssociationRegistry(
            global_constraints=self.global_constraints
        )
        pool = AssociationPool.read(self.pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            asn_name = asn.asn_name
            m = re.match(LEVEL3_ASN_NAME_REGEX, asn_name)
            yield helpers.not_none, m
            yield helpers.check_equal, m.groupdict()['acid'], '0002'
