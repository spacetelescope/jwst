import pytest
import re

from .helpers import (
    combine_pools,
    func_fixture,
    generate_params,
    registry_level3_only,
    t_path,
)

from .. import (AssociationPool, generate)
from ..lib.dms_base import DMSAttrConstraint

# Temporarily skip if running under Travis
# pytestmark = pytest.mark.skipif(
#     "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
#     reason='Temporarily disable due to performance issues'
# )

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
EMPTY = (None, '', 'NULL', 'Null', 'null', 'CLEAR', 'Clear', 'clear', 'F', 'f', 'N', 'n')


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


def test_level35_names(pool_file):
    rules = registry_level3_only()
    pool = AssociationPool.read(pool_file)
    asns = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] == 'Asn_IFU':
            m = re.match(LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX, product_name)
        else:
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None


def test_level3_names(pool_file, global_constraints):
    rules = registry_level3_only(
        global_constraints=global_constraints
    )
    pool = AssociationPool.read(pool_file)
    asns = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] == 'Asn_IFU':
            m = re.match(LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX, product_name)
        else:
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None
        assert m.groupdict()['acid'] == 'o002'


@pytest.mark.xfail(
    reason='Unknown, need to investigate',
    run=False
)
def test_multiple_optelems(pool_file):
    rules = registry_level3_only()
    pool = AssociationPool.read(pool_file)
    asns = generate(pool, rules)
    for asn in asns:
        product_name = asn['products'][0]['name']
        if asn['asn_rule'] != 'Asn_IFU':
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
            assert m is not None
            try:
                value = '-'.join(asn.constraints['opt_elem2']['found_values'])
            except KeyError:
                value = None
            if value in EMPTY:
                assert '-' not in m.groupdict()['opt_elem']
            else:
                assert '-' in m.groupdict()['opt_elem']
