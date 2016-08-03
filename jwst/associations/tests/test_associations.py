"""test_associations: Test of general Association functionality."""
from __future__ import absolute_import

import nose.tools as nt

from . import helpers
from .helpers import full_pool_rules

from .. import (
    Association,
    AssociationError,
    AssociationRegistry,
    generate)
from ..registry import (
    import_from_file,
    find_member
)


class TestAssociations():

    pools_size = [
        (
            helpers.t_path('data/jw93060_20150312T160130_pool.csv'),
            14
        ),
    ]

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # Basic Association object
    def test_read_assoc_defs(self):
        rules = AssociationRegistry(
            [helpers.t_path('data/asn_rules_set1.py')],
            include_default=False
        )
        assert len(rules) >= 2
        assert 'DMS_Level3_Base_Set1' not in rules
        valid_rules = ['Asn_Dither_Set1', 'Asn_WFS_Set1']
        for rule in valid_rules:
            yield helpers.check_in_list, rule, rules

    def test_read_assoc_defs_fromdefault(self):
        rules = AssociationRegistry()
        assert len(rules) >= 3
        assert 'DMS_Level3_Base' not in rules
        valid_rules = ['Asn_Image', 'Asn_WFSCMB']
        for rule in valid_rules:
            yield helpers.check_in_list, rule, rules

    @nt.raises(AssociationError)
    def test_nodefs(self):
        rules = AssociationRegistry(include_default=False)

    def test_multi_rules(self):
        rule_files = [
            helpers.t_path('data/asn_rules_set1.py'),
            helpers.t_path('data/asn_rules_set2.py')
        ]
        rules = AssociationRegistry(rule_files, include_default=False)
        assert len(rules) == 4
        assert 'DMS_Level3_Base_Set1' not in rules
        assert 'DMS_Level3_Base_Set2' not in rules
        valid_rules = [
            'Asn_Dither_Set1',
            'Asn_Dither_Set2',
            'Asn_WFS_Set1',
            'Asn_WFS_Set2'
        ]

        for rule in valid_rules:
            yield helpers.check_in_list, rule, rules

    def test_base_instatiation(self):
        """Create an association without any initialization"""
        assert Association()

    def test_global_constraints(self, full_pool_rules):
        """Test that global constraints get applied to all rules"""
        full_pool, default_rules = full_pool_rules

        tests = {
            'exists': {
                'constraints': {
                    'obs_id': {
                        'value': 'V99009001001P0000000002101',
                        'inputs': ['OBS_ID']
                    }
                },
                'pool': full_pool,
                'n_asns': 1,
            },
            'empty': {
                'constraints': {
                    'obs_id': {
                        'value': 'junk',
                        'inputs': ['OBS_ID']
                    }
                },
                'pool': helpers.t_path('data/pool_candidates.csv'),
                'n_asns': 0,
            },
            'combined_candidates': {
                'constraints': {
                    'asn_candidate_id': {
                        'value': '.+(o001|o002).+',
                        'inputs': ['ASN_CANDIDATE'],
                        'force_unique': False,
                    }
                },
                'pool': helpers.t_path('data/pool_candidates.csv'),
                'n_asns': 2,
            },
            'exclusive_candidates': {
                'constraints': {
                    'asn_candidate_id': {
                        'value': '.+(o001|o002).+',
                        'inputs': ['ASN_CANDIDATE'],
                        'force_unique': True,
                    }
                },
                'pool': helpers.t_path('data/pool_candidates.csv'),
                'n_asns': 4,
            },
        }

        for test_name, test in tests.items():
            rules = AssociationRegistry(
                global_constraints=test['constraints']
            )
            assert len(rules) >= 3
            for constraint in test['constraints']:
                for rule in rules:
                    assert constraint in rules[rule].GLOBAL_CONSTRAINTS

            pool = helpers.combine_pools(test['pool'])
            asns, orphaned = generate(pool, rules)
            assert len(asns) == test['n_asns']

    def test_rulesets(self):
        """Test finding members in a ruleset"""
        rule_file = helpers.t_path('../lib/association_rules.py')
        module = import_from_file(rule_file)
        schemas = [schema for schema in find_member(module, 'ASN_SCHEMA')]
        assert len(schemas) == 2
