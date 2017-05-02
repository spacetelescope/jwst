"""test_associations: Test of general Association functionality."""
from __future__ import absolute_import
import os
import pytest

from . import helpers

from .. import (
    Association,
    AssociationError,
    AssociationRegistry,
    generate)
from ..registry import (
    import_from_file,
    find_member
)

# Temporarily skip if running under Travis
pytestmark = pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason='Temporarily disable due to performance issues'
)


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
    valid_rules = ['Asn_Image', 'Asn_WFSCMB']
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


def test_base_instatiation():
    """Create an association without any initialization"""
    assert Association()


@pytest.mark.parametrize(
    'constraints, pool, n_asns',
    [
        (
            {
                'obs_id': {
                    'value': 'V99009001001P0000000002101',
                    'inputs': ['OBS_ID']
                }
            },
            helpers.t_path('data/mega_pool.csv'),
            2,
        ),
        (
            {
                'obs_id': {
                    'value': 'junk',
                    'inputs': ['OBS_ID']
                }
            },
            helpers.t_path('data/pool_001_candidates.csv'),
            0,
        ),
        (
            {
                'asn_candidate_id': {
                    'value': '.+(o001|o002).+',
                    'inputs': ['ASN_CANDIDATE'],
                    'force_unique': False,
                }
            },
            helpers.t_path('data/pool_001_candidates.csv'),
            3,
        ),
        (
            {
                'asn_candidate_id': {
                    'value': '.+(o001|o002).+',
                    'inputs': ['ASN_CANDIDATE'],
                    'force_unique': True,
                }
            },
            helpers.t_path('data/pool_001_candidates.csv'),
            5,
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
            assert constraint in rules[rule].GLOBAL_CONSTRAINTS

    pool = helpers.combine_pools(pool)
    asns, orphaned = generate(pool, rules)
    assert len(asns) == n_asns


def test_rulesets():
    """Test finding members in a ruleset"""
    rule_file = helpers.t_path('../lib/association_rules.py')
    module = import_from_file(rule_file)
    schemas = [schema for schema in find_member(module, 'ASN_SCHEMA')]
    assert len(schemas) == 2
