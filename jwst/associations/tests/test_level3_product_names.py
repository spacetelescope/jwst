"""Check formatting of the level 3 product names"""

import re

import pytest
from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests.helpers import (
    combine_pools,
    registry_level3_only,
)
from jwst.associations import generate
from jwst.associations.lib.dms_base import DMSAttrConstraint


LEVEL3_PRODUCT_NAME_REGEX = (
    r"jw"
    r"(?P<program>\d{5})"
    r"-(?P<acid>[a-z]\d{3,4})"
    r"_(?P<target>(?:t\d{3})|(?:\{source_id\}))"
    r"(?:-(?P<epoch>epoch\d+))?"
    r"_(?P<instrument>.+?)"
    r"_(?P<opt_elem>.+)"
)

LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX = (
    r"jw"
    r"(?P<program>\d{5})"
    r"-(?P<acid>[a-z]\d{3,4})"
    r"_(?P<target>(?:t\d{3})|(?:s\d{5}))"
    r"(?:-(?P<epoch>epoch\d+))?"
    r"_(?P<instrument>.+?)"
)

LEVEL3_PRODUCT_NAME_NRS_FSS_REGEX = (
    r"jw"
    r"(?P<program>\d{5})"
    r"-(?P<acid>[a-z]\d{3,4})"
    r"_(?P<target>(?:t\d{3})"
    r"-(?:\{source_id\}))"
    r"(?:-(?P<epoch>epoch\d+))?"
    r"_(?P<instrument>.+?)"
    r"_(?P<opt_elem>.+)"
)

# Null values
EMPTY = (None, "", "NULL", "Null", "null", "F", "f", "N", "n")


@pytest.fixture(scope="module")
def pool_file():
    return get_pkg_data_filename(
        "data/pool_018_all_exptypes.csv", package="jwst.associations.tests"
    )


@pytest.fixture(scope="module")
def global_constraints():
    constraint = DMSAttrConstraint(
        name="asn_candidate",
        value=[".+o002.+"],
        sources=["asn_candidate"],
        force_unique=True,
        is_acid=True,
        evaluate=True,
    )
    return constraint


def test_level3_productname_components_discovered():
    rules = registry_level3_only()
    pool = combine_pools(
        get_pkg_data_filename("data/pool_002_image_miri.csv", package="jwst.associations.tests")
    )
    asn = generate(pool, rules)[0]
    match = re.match(LEVEL3_PRODUCT_NAME_REGEX, asn["products"][0]["name"])
    assert match is not None
    matches = match.groupdict()
    assert matches["program"] == "99009"
    assert matches["acid"] == "a3001"
    assert matches["target"] == "t001"
    assert matches["instrument"] == "miri"
    assert matches["opt_elem"] == "f560w"


def test_level3_productname_components_acid():
    global_constraints = DMSAttrConstraint(
        name="asn_candidate_ids",
        value=".+o001.+",
        sources=["asn_candidate"],
        force_unique=True,
        is_acid=True,
        evaluate=True,
    )
    rules = registry_level3_only(global_constraints=global_constraints)
    pool = combine_pools(
        get_pkg_data_filename("data/pool_002_image_miri.csv", package="jwst.associations.tests")
    )
    asn = generate(pool, rules)[0]
    match = re.match(LEVEL3_PRODUCT_NAME_REGEX, asn["products"][0]["name"])
    assert match is not None
    matches = match.groupdict()
    assert matches["program"] == "99009"
    assert matches["acid"] == "o001"
    assert matches["target"] == "t001"
    assert matches["instrument"] == "miri"
    assert matches["opt_elem"] == "f560w"


def test_level3_names(pool_file, global_constraints):
    rules = registry_level3_only(global_constraints=global_constraints)
    pool = combine_pools(pool_file)
    asns = generate(pool, rules)
    for asn in asns:
        product_name = asn["products"][0]["name"]
        if asn["asn_rule"] == "Asn_Lv3MIRMRS":
            m = re.match(LEVEL3_PRODUCT_NAME_NO_OPTELEM_REGEX, product_name)
        else:
            m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None
        assert m.groupdict()["acid"] == "o002"


def test_multiple_optelems(pool_file):
    rules = registry_level3_only()
    pool = combine_pools(pool_file)
    asns = generate(pool, rules, finalize=False)
    for asn in asns:
        product_name = asn["products"][0]["name"]
        if asn["asn_rule"] == "Asn_Lv3MIRMRS":
            continue
        if asn["asn_rule"] == "Asn_Lv3NRSFSS":
            regex_to_match = LEVEL3_PRODUCT_NAME_NRS_FSS_REGEX
        else:
            regex_to_match = LEVEL3_PRODUCT_NAME_REGEX

        m = re.match(regex_to_match, product_name)
        assert m is not None

        # there should always be an opt_elem
        values = ["-".join(asn.constraints["opt_elem"].found_values)]

        # there may also be an opt_elem2, fixed slit or 2, or a subarray
        for extra in ["opt_elem2", "fxd_slit", "fxd_slit2", "subarray"]:
            # special rules for fixed slit for NRS FS:
            # it gets a format placeholder instead of the value
            if asn["asn_rule"] == "Asn_Lv3NRSFSS":
                if extra == "fxd_slit":
                    values.append("{slit_name}")
                    continue
                elif extra == "fxd_slit2":
                    continue

            try:
                value = "-".join(asn.constraints[extra].found_values)
            except KeyError:
                value = None

            # empty values and subarray = full are not recorded
            if value not in EMPTY and value != "full":
                values.append(value)

        assert m.groupdict()["opt_elem"] == "-".join(values)


def test_tso3_names():
    rules = registry_level3_only()
    tso_pool = get_pkg_data_filename("data/pool_021_tso.csv", package="jwst.associations.tests")
    pool = combine_pools(tso_pool)
    asns = generate(pool, rules, finalize=False)
    for asn in asns:
        product_name = asn["products"][0]["name"]

        m = re.match(LEVEL3_PRODUCT_NAME_REGEX, product_name)
        assert m is not None

        # there should always be an opt_elem
        values = ["-".join(asn.constraints["opt_elem"].found_values)]

        # there may also be an opt_elem2, fixed slit or 2, or a subarray
        for extra in ["opt_elem2", "fxd_slit", "fxd_slit2", "subarray"]:
            try:
                value = "-".join(asn.constraints[extra].found_values)
            except KeyError:
                value = None

            # empty values and subarray = full are not recorded
            if value not in EMPTY and value != "full":
                values.append(value)

        assert m.groupdict()["opt_elem"] == "-".join(values)
