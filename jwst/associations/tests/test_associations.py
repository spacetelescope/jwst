"""test_associations: Test of general Association functionality."""

import pytest
from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests import helpers
from jwst.associations import Association, AssociationError, AssociationRegistry, generate
from jwst.associations.lib.dms_base import DMSAttrConstraint


# Basic Association object
def test_read_assoc_defs():
    rules = AssociationRegistry(
        [get_pkg_data_filename("data/asn_rules_set1.py", package="jwst.associations.tests")],
        include_default=False,
    )
    assert len(rules) >= 2
    rule_names = helpers.get_rule_names(rules)
    assert "DMSLevel3BaseSet1" not in rules
    valid_rules = ["AsnDitherSet1", "AsnWFSSet1"]
    for rule in valid_rules:
        assert rule in rule_names


def test_read_assoc_defs_fromdefault():
    rules = AssociationRegistry()
    assert len(rules) >= 3
    rule_names = helpers.get_rule_names(rules)
    assert "DMS_Level3_Base" not in rules
    valid_rules = ["Asn_Lv3Image", "Asn_Lv3WFSCMB"]
    for rule in valid_rules:
        assert rule in rule_names


def test_registry_backref():
    rules = AssociationRegistry()
    for name, rule in rules.items():
        assert rule.registry == rules


def test_nodefs():
    with pytest.raises(AssociationError):
        AssociationRegistry(include_default=False)


def test_multi_rules():
    rule_files = [
        get_pkg_data_filename("data/asn_rules_set1.py", package="jwst.associations.tests"),
        get_pkg_data_filename("data/asn_rules_set2.py", package="jwst.associations.tests"),
    ]
    rules = AssociationRegistry(rule_files, include_default=False)
    assert len(rules) == 4
    rule_names = helpers.get_rule_names(rules)
    assert "DMSLevel3BaseSet1" not in rule_names
    assert "DMSLevel3BaseSet2" not in rule_names
    valid_rules = ["AsnDitherSet1", "AsnDitherSet2", "AsnWFSSet1", "AsnWFSSet2"]

    for rule in valid_rules:
        assert rule in rule_names


def test_base_instantiation():
    """Create an association without any initialization"""
    assert Association()


@pytest.mark.parametrize(
    "constraints, pool, n_asns",
    [
        (
            DMSAttrConstraint(
                name="obs_id", value="V99009001001P0000000002101", sources=["obs_id"]
            ),
            get_pkg_data_filename(
                "data/pool_018_all_exptypes.csv", package="jwst.associations.tests"
            ),
            1,
        ),
        (
            DMSAttrConstraint(name="obs_id", value="junk", sources=["obs_id"]),
            get_pkg_data_filename(
                "data/pool_001_candidates.csv", package="jwst.associations.tests"
            ),
            0,
        ),
        (
            DMSAttrConstraint(
                name="asn_candidate_id",
                value=".+(o001|o002).+",
                sources=["asn_candidate"],
                force_unique=False,
                is_acid=True,
                evaluate=True,
            ),
            get_pkg_data_filename(
                "data/pool_001_candidates.csv", package="jwst.associations.tests"
            ),
            22,
        ),
        (
            DMSAttrConstraint(
                name="asn_candidate_id",
                value=".+(o001|o002).+",
                sources=["asn_candidate"],
                force_unique=True,
                is_acid=True,
                evaluate=True,
            ),
            get_pkg_data_filename(
                "data/pool_001_candidates.csv", package="jwst.associations.tests"
            ),
            24,
        ),
    ],
)
def test_global_constraints(constraints, pool, n_asns):
    """Test that global constraints get applied to all rules"""
    rules = AssociationRegistry(global_constraints=constraints)
    assert len(rules) >= 3
    for constraint in constraints:
        for rule in rules:
            assert constraint in rules[rule].GLOBAL_CONSTRAINT

    pool = helpers.combine_pools(pool)
    asns = generate(pool, rules)
    assert len(asns) == n_asns
