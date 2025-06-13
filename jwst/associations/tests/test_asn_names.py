import pytest
import re

from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests import helpers
from jwst.associations import generate
from jwst.associations import load_asn
from jwst.associations.lib.utilities import constrain_on_candidates

LEVEL3_ASN_ACID_NAME_REGEX = (
    r"jw"
    r"(?P<program>\d{5})"
    r"-(?P<acid>(o|c)\d{3,4})"
    r"_(?P<asn_type>\w+)"
    r"_(?P<sequence>\d{5})"
    r"_asn"
)
LEVEL3_ASN_DISCOVERED_NAME_REGEX = (
    r"jw"
    r"(?P<program>\d{5})"
    r"-(?P<acid>a\d{4})"
    r"_(?P<asn_type>\w+)"
    r"_(?P<sequence>\d{5})"
    r"_asn"
)

LEVEL3_ASN_WITH_VERSION = (
    r"jw"
    r"(?P<program>\d{5})"
    r"-(?P<acid>[a-z]\d{3,4})"
    r"_(?P<stamp>.+)"
    r"_(?P<asn_type>.+)"
    r"_(?P<sequence>\d{5})"
    r"_asn"
)

all_candidates = constrain_on_candidates(None)


def test_level3_asn_names():
    pool_path = get_pkg_data_filename(
        "data/pool_020_00009_image_miri.csv", package="jwst.associations.tests"
    )
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level3_only()
    asns = generate(pool, rules)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        if any(getattr(c, "is_acid", False) for c in asn.constraints):
            m = re.match(LEVEL3_ASN_ACID_NAME_REGEX, name)
        else:
            m = re.match(LEVEL3_ASN_DISCOVERED_NAME_REGEX, name)
        assert m is not None


def test_level3_asn_names_with_version():
    pool_path = get_pkg_data_filename(
        "data/pool_020_00009_image_miri.csv", package="jwst.associations.tests"
    )
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level3_only()
    asns = generate(pool, rules, version_id=True)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        m = re.match(LEVEL3_ASN_WITH_VERSION, name)
        assert m is not None


def test_level2_asn_names():
    pool_path = get_pkg_data_filename(
        "data/pool_020_00009_image_miri.csv", package="jwst.associations.tests"
    )
    pool = helpers.combine_pools(pool_path)
    rules = helpers.registry_level2_only(global_constraints=all_candidates)
    asns = generate(pool, rules)
    assert len(asns) > 0
    for asn in asns:
        name = asn.asn_name
        if any(getattr(c, "is_acid", False) for c in asn.constraints):
            m = re.match(LEVEL3_ASN_ACID_NAME_REGEX, name)
        else:
            m = re.match(LEVEL3_ASN_DISCOVERED_NAME_REGEX, name)
        assert m is not None


def test_level2_asn_names_with_version():
    pool_path = get_pkg_data_filename(
        "data/pool_020_00009_image_miri.csv", package="jwst.associations.tests"
    )
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
    fname = get_pkg_data_filename(
        "data/asn_level2_bad_path.json", package="jwst.associations.tests"
    )
    with open(fname) as fd:
        with pytest.warns(UserWarning):
            load_asn(fd)
