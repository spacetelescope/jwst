import pytest
import re

from jwst.associations.tests.helpers import (
    func_fixture,
    generate_params,
    registry_level3_only,
    t_path,
)

from jwst.associations import (AssociationPool, generate)
from jwst.associations.lib.dms_base import DMSAttrConstraint


LEVEL3_PRODUCT_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>[a-z]\d{3,4})'
    '_(?P<target>(?:t\d{3})|(?:\{source_id\}))'
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
EMPTY = (None, '', 'NULL', 'Null', 'null', 'F', 'f', 'N', 'n')


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


@pytest.mark.slow
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
