from __future__ import absolute_import

import os
import re

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)

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

    pool_file = helpers.t_path(
        'data/pool_001_candidates.csv'
    )

    global_constraints = {
        'asn_candidate_ids': {
            'value': ['o002'],
            'inputs': ['ASN_CANDIDATE'],
            'force_unique': True,
        }
    }

    def test_level35_names(self):
        rules = AssociationRegistry()
        pool = helpers.combine_pools(self.pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            asn_name = asn.asn_name
            m = re.match(LEVEL3_ASN_NAME_REGEX, asn_name)
            yield helpers.not_none, m

    def test_level3_names(self):
        rules = AssociationRegistry(
            global_constraints=self.global_constraints
        )
        pool = helpers.combine_pools(self.pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            asn_name = asn.asn_name
            m = re.match(LEVEL3_ASN_NAME_REGEX, asn_name)
            yield helpers.not_none, m
            yield helpers.check_equal, m.groupdict()['acid'], 'o002'
