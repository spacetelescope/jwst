from __future__ import absolute_import

import pytest
import re

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)

LEVEL3_PRODUCT_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>[a-z]\d{3,4})'
    '_(?P<target>(?:t\d{3})|(?:s\d{5}))'
    '(?:-(?P<epoch>epoch\d+))?'
    '_(?P<instrument>.+?)'
    '_(?P<opt_elem>.+)'
)


# Null values
EMPTY = (None, 'NULL', 'CLEAR')


pool_file = helpers.func_fixture(
    helpers.generate_params,
    scope='module',
    params=[
        helpers.t_path('data/mega_pool.csv'),
    ]
)


global_constraints = helpers.func_fixture(
    helpers.generate_params,
    scope='module',
    params=[
        {
            'asn_candidate': {
                'value': ['.+o002.+'],
                'inputs': ['ASN_CANDIDATE'],
                'force_unique': True,
                'is_acid': True,
            }
        },
    ]
)


class TestProductNames():

    def test_level35_names(self, pool_file):
        rules = AssociationRegistry()
        pool = AssociationPool.read(pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            product_name = asn['products'][0]['name']
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            assert m is not None

    def test_level3_names(self, pool_file, global_constraints):
        rules = AssociationRegistry(
            global_constraints=global_constraints
        )
        pool = AssociationPool.read(pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            product_name = asn['products'][0]['name']
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            assert m is not None
            assert m.groupdict()['acid'] == 'o002'

    def test_multiple_optelems(self, pool_file):
        rules = AssociationRegistry()
        pool = AssociationPool.read(pool_file)
        (asns, orphaned) = generate(pool, rules)
        for asn in asns:
            product_name = asn['products'][0]['name']
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            assert m is not None
            if asn.constraints['opt_elem2']['value'] in EMPTY:
                assert '-' not in m.groupdict()['opt_elem']
            else:
                assert '-' in m.groupdict()['opt_elem']
