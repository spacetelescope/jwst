"""Check formatting of the level 3 product names"""
import re

import pytest

from jwst.associations.tests.helpers import (
    combine_pools,
    registry_level3_only,
    t_path,
)

from jwst.associations import (AssociationPool, generate)
from jwst.associations.lib.dms_base import DMSAttrConstraint


LEVEL3_PRODUCT_NAME_REGEX = (
    r'jw'
    r'(?P<program>\d{5})'
    r'-(?P<acid>[a-z]\d{3,4})'
    r'_(?P<target>(?:t\d{3})|(?:\{source_id\}))'
    r'(?:-(?P<epoch>epoch\d+))?'
    r'_(?P<instrument>.+?)'
    r'_(?P<opt_elem>.+)'
)

LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX = (
    r'jw'
    r'(?P<program>\d{5})'
    r'-(?P<acid>[a-z]\d{3,4})'
    r'_(?P<target>(?:t\d{3})|(?:s\d{5}))'
    r'(?:-(?P<epoch>epoch\d+))?'
    r'_(?P<instrument>.+?)'
)

# Null values
EMPTY = (None, '', 'NULL', 'Null', 'null', 'F', 'f', 'N', 'n')


@pytest.fixture(scope='module')
def pool_file():
    return t_path('data/pool_018_all_exptypes.csv')


@pytest.fixture(scope='module')
def global_constraints():
    constraint = DMSAttrConstraint(
        name='asn_candidate',
        value=['.+o002.+'],
        sources=['asn_candidate'],
        force_unique=True,
        is_acid=True,
        evaluate=True,
    )
    return constraint


def test_level3_productname_components_discovered():
    rules = registry_level3_only()
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns = generate(pool, rules)
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
    global_constraints = DMSAttrConstraint(
        name='asn_candidate_ids',
        value='.+o001.+',
        sources=['asn_candidate'],
        force_unique=True,
        is_acid=True,
        evaluate=True,
    )
    rules = registry_level3_only(global_constraints=global_constraints)
    pool = combine_pools(t_path('data/pool_002_image_miri.csv'))
    asns = generate(pool, rules)
    asn = asns[0]
    match = re.match(LEVEL3_PRODUCT_NAME_REGEX, asn['products'][0]['name'])
    assert match is not None
    matches = match.groupdict()
    assert matches['program'] == '99009'
    assert matches['acid'] == 'o001'
    assert matches['target'] == 't001'
    assert matches['instrument'] == 'miri'
    assert matches['opt_elem'] == 'f560w'


def test_level3_names(pool_file, global_constraints):
    rules = registry_level3_only(
        global_constraints=global_constraints
    )
    pool = AssociationPool.read(pool_file)
    asns = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] == 'Asn_Lv3MIRMRS':
            m = re.match(LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX, product_name)
        else:
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None
        assert m.groupdict()['acid'] == 'o002'


def test_multiple_optelems(pool_file):
    rules = registry_level3_only()
    pool = AssociationPool.read(pool_file)
    asns = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] != 'Asn_Lv3MIRMRS':
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            assert m is not None
            try:
                value = '-'.join(asn.constraints['opt_elem2'].found_values)
            except KeyError:
                value = None
            if value in EMPTY:
                assert '-' not in m.groupdict()['opt_elem']
            else:
                assert '-' in m.groupdict()['opt_elem']
