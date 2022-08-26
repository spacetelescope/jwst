import os
import pytest
import re
import warnings

from jwst.associations.tests import helpers

from jwst import associations
from jwst.associations import generate
from jwst.associations import load_asn
from jwst.associations.main import constrain_on_candidates

LEVEL3_ASN_ACID_NAME_REGEX = (
    r'jw'
    r'(?P<program>\d{5})'
    r'-(?P<acid>(o|c)\d{3,4})'
    r'_(?P<asn_type>\w+)'
    r'_(?P<sequence>\d{3})'
    r'_asn'
)
LEVEL3_ASN_DISCOVERED_NAME_REGEX = (
    r'jw'
    r'(?P<program>\d{5})'
    r'-(?P<acid>a\d{4})'
    r'_(?P<asn_type>\w+)'
    r'_(?P<sequence>\d{3})'
    r'_asn'
)

LEVEL3_ASN_WITH_VERSION = (
    r'jw'
    r'(?P<program>\d{5})'
    r'-(?P<acid>[a-z]\d{3,4})'
    r'_(?P<stamp>.+)'
    r'_(?P<asn_type>.+)'
    r'_(?P<sequence>\d{3})'
    r'_asn'
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


def test_bad_expnames():
    """
    Ensure warning gets raised during load_asn when the association file
    contains path data in the expname.
    """
    bad_asn = "./data/asn_level2_bad_path.json"
    fname = os.path.abspath(os.path.join(os.path.dirname(__file__), bad_asn))
    with open(fname) as fd:
        with pytest.warns(UserWarning):
            asn = load_asn(fd)
