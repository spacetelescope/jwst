"""test_associations: Test of general Association functionality."""
import pytest

from jwst.associations.tests import helpers

from jwst.associations import (
    Association,
    AssociationError,
    AssociationRegistry,
    generate)
from jwst.associations.registry import (
    import_from_file,
)
from jwst.associations.lib.dms_base import DMSAttrConstraint


# Basic Association object
def test_read_assoc_defs():
    rules = AssociationRegistry(
        [helpers.t_path('data/asn_rules_set1.py')],
        include_default=False
    )
    assert len(rules) >= 2
    rule_names = helpers.get_rule_names(rules)
    assert 'DMS_Level3_Base_Set1' not in rules
    valid_rules = ['Asn_Dither_Set1', 'Asn_WFS_Set1']
    for rule in valid_rules:
        assert rule in rule_names


def test_read_assoc_defs_fromdefault():
    rules = AssociationRegistry()
    assert len(rules) >= 3
    rule_names = helpers.get_rule_names(rules)
    assert 'DMS_Level3_Base' not in rules
    valid_rules = ['Asn_Lv3Image', 'Asn_Lv3WFSCMB']
    for rule in valid_rules:
        assert rule in rule_names


def test_registry_backref():
    rules = AssociationRegistry()
    for name, rule in rules.items():
        assert rule.registry == rules


def test_nodefs():
    with pytest.raises(AssociationError):
        rules = AssociationRegistry(include_default=False)


def test_multi_rules():
    rule_files = [
        helpers.t_path('data/asn_rules_set1.py'),
        helpers.t_path('data/asn_rules_set2.py')
    ]
    rules = AssociationRegistry(rule_files, include_default=False)
    assert len(rules) == 4
    rule_names = helpers.get_rule_names(rules)
    assert 'DMS_Level3_Base_Set1' not in rule_names
    assert 'DMS_Level3_Base_Set2' not in rule_names
    valid_rules = [
        'Asn_Dither_Set1',
        'Asn_Dither_Set2',
        'Asn_WFS_Set1',
        'Asn_WFS_Set2'
    ]

    for rule in valid_rules:
        assert rule in rule_names


def test_base_instantiation():
    """Create an association without any initialization"""
    assert Association()


@pytest.mark.parametrize(
    'constraints, pool, n_asns',
    [
        (
            DMSAttrConstraint(
                name='obs_id',
                value='V99009001001P0000000002101',
                sources=['obs_id']
            ),
            helpers.t_path('data/pool_018_all_exptypes.csv'),
            1,
        ),
        (
            DMSAttrConstraint(
                name='obs_id',
                value='junk',
                sources=['obs_id']
            ),
            helpers.t_path('data/pool_001_candidates.csv'),
            0,
        ),
        (
            DMSAttrConstraint(
                name='asn_candidate_id',
                value='.+(o001|o002).+',
                sources=['asn_candidate'],
                force_unique=False,
                is_acid=True,
                evaluate=True,
            ),
            helpers.t_path('data/pool_001_candidates.csv'),
            22,
        ),
        (
            DMSAttrConstraint(
                name='asn_candidate_id',
                value='.+(o001|o002).+',
                sources=['asn_candidate'],
                force_unique=True,
                is_acid=True,
                evaluate=True,
            ),
            helpers.t_path('data/pool_001_candidates.csv'),
            24,
        ),
    ]
)
def test_global_constraints(constraints, pool, n_asns):
    """Test that global constraints get applied to all rules"""
    rules = AssociationRegistry(
        global_constraints=constraints
    )
    assert len(rules) >= 3
    for constraint in constraints:
        for rule in rules:
            assert constraint in rules[rule].GLOBAL_CONSTRAINT

    pool = helpers.combine_pools(pool)
    asns = generate(pool, rules)
    assert len(asns) == n_asns
