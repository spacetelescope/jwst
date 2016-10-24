from __future__ import absolute_import

import re

from .helpers import (
    combine_pools,
    func_fixture,
    generate_params,
    t_path,
)

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

LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>[a-z]\d{3,4})'
    '_(?P<target>(?:t\d{3})|(?:s\d{5}))'
    '(?:-(?P<epoch>epoch\d+))?'
    '_(?P<instrument>.+?)'
)

# Null values
EMPTY = (None, 'NULL', 'CLEAR')


pool_file = func_fixture(
    generate_params,
    scope='module',
    params=[
        t_path('data/mega_pool.csv'),
    ]
)


global_constraints = func_fixture(
    generate_params,
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


def test_level3_productname_components_discovered():
    rules = AssociationRegistry()
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns, orphaned = generate(pool, rules)
    asn = asns[0]
    match = re.match(LEVEL3_PRODUCT_NAME_REGEX, asn['products'][0]['name'])
    assert match is not None
    matches = match.groupdict()
    assert matches['program'] == '99009'
    assert matches['acid'] == 'a3001'
    assert matches['target'] == 't001'
    assert matches['instrument'] == 'miri'
    assert matches['opt_elem'] == 'f560w'


def test_level3_productname_components_acid():
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
    match = re.match(LEVEL3_PRODUCT_NAME_REGEX, asn['products'][0]['name'])
    assert match is not None
    matches = match.groupdict()
    assert matches['program'] == '99009'
    assert matches['acid'] == 'o001'
    assert matches['target'] == 't001'
    assert matches['instrument'] == 'miri'
    assert matches['opt_elem'] == 'f560w'


def test_level35_names(pool_file):
    rules = AssociationRegistry()
    pool = AssociationPool.read(pool_file)
    (asns, orphaned) = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] == 'Asn_MIRI_IFU':
            m = re.match(LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX, product_name)
        else:
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None


def test_level3_names(pool_file, global_constraints):
    rules = AssociationRegistry(
        global_constraints=global_constraints
    )
    pool = AssociationPool.read(pool_file)
    (asns, orphaned) = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] == 'Asn_MIRI_IFU':
            m = re.match(LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX, product_name)
        else:
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None
        assert m.groupdict()['acid'] == 'o002'


def test_multiple_optelems(pool_file):
    rules = AssociationRegistry()
    pool = AssociationPool.read(pool_file)
    (asns, orphaned) = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] != 'Asn_MIRI_IFU':
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            assert m is not None
            try:
                value = asn.constraints['opt_elem2']['value']
            except KeyError:
                value = None
            if value in EMPTY:
                assert '-' not in m.groupdict()['opt_elem']
            else:
                assert '-' in m.groupdict()['opt_elem']
