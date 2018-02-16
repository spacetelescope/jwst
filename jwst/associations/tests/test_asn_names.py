import pytest
import re

from . import helpers

from .. import generate
from ..main import constrain_on_candidates

LEVEL3_ASN_ACID_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>(o|c)\d{3,4})'
    '_(?P<asn_type>\w+)'
    '_(?P<sequence>\d{3})'
    '_asn'
)
LEVEL3_ASN_DISCOVERED_NAME_REGEX = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>a\d{4})'
    '_(?P<asn_type>\w+)'
    '_(?P<sequence>\d{3})'
    '_asn'
)

LEVEL3_ASN_WITH_VERSION = (
    'jw'
    '(?P<program>\d{5})'
    '-(?P<acid>[a-z]\d{3,4})'
    '_(?P<stamp>.+)'
    '_(?P<asn_type>.+)'
    '_(?P<sequence>\d{3})'
    '_asn'
)

all_candidates = constrain_on_candidates(None)

pool_params = pytest.fixture(
    scope='module',
    params=[
        'data/pool_020_00009_image_miri.csv'
    ]
)(helpers.generate_params)


def test_level3_asn_names(pool_params):
    pool_path = helpers.t_path(pool_params)
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level3_only()
    asns = generate(pool, rules)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        if any(
                getattr(c, 'is_acid', False)
                for c in asn.constraints
        ):
            m = re.match(LEVEL3_ASN_ACID_NAME_REGEX, name)
        else:
            m = re.match(LEVEL3_ASN_DISCOVERED_NAME_REGEX, name)
        assert m is not None


def test_level3_asn_names_with_version(pool_params):
    pool_path = helpers.t_path(pool_params)
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level3_only()
    asns = generate(pool, rules, version_id=True)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        m = re.match(LEVEL3_ASN_WITH_VERSION, name)
        assert m is not None


def test_level2_asn_names(pool_params):
    pool_path = helpers.t_path(pool_params)
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level2_only(global_constraints=all_candidates)
    asns = generate(pool, rules)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        if any(
                getattr(c, 'is_acid', False)
                for c in asn.constraints
        ):
            m = re.match(LEVEL3_ASN_ACID_NAME_REGEX, name)
        else:
            m = re.match(LEVEL3_ASN_DISCOVERED_NAME_REGEX, name)
        assert m is not None


def test_level2_asn_names_with_version(pool_params):
    pool_path = helpers.t_path(pool_params)
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level2_only(global_constraints=all_candidates)
    asns = generate(pool, rules, version_id=True)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        m = re.match(LEVEL3_ASN_WITH_VERSION, name)
        assert m is not None
