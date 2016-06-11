"""test_associations: Test of general Association functionality."""
from __future__ import absolute_import

import nose.tools as nt
from nose import SkipTest

from . import helpers

from .. import (
    AssociationError,
    AssociationRegistry,
    AssociationPool,
    generate)


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
            [helpers.t_path('tests/data/asn_rules_set1.py')],
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
        valid_rules = ['Asn_Dither', 'Asn_WFSCMB']
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
        raise SkipTest('Need to implement.')

    def test_global_constraints(self):
        """Test that global constraints get applied to all rules"""

        tests = {
            'exists': {
                'constraints': {
                    'obs_id': {
                        'value': 'V93060001004P0000000001101',
                        'inputs': ['OBS_ID']
                    }
                },
                'pool': helpers.t_path('data/jw93060_20150312T160130_pool.csv'),
                'n_asns': 6,
            },
            'empty': {
                'constraints': {
                    'obs_id': {
                        'value': 'junk',
                        'inputs': ['OBS_ID']
                    }
                },
                'pool': helpers.t_path('data/jw93060_20150312T160130_pool.csv'),
                'n_asns': 0,
            },
            'combined_candidates': {
                'constraints': {
                    'asn_candidate_id': {
                        'value': '1|2',
                        'inputs': ['ASN_CANDIDATE_ID', 'OBS_NUM'],
                        'force_unique': False,
                    }
                },
                'pool': helpers.t_path('data/jw93060_002_20150312T160130_pool.csv'),
                'n_asns': 6,
            },
            'exclusive_candidates': {
                'constraints': {
                    'asn_candidate_id': {
                        'value': '1|2',
                        'inputs': ['ASN_CANDIDATE_ID', 'OBS_NUM'],
                        'force_unique': True,
                    }
                },
                'pool': helpers.t_path('data/jw93060_002_20150312T160130_pool.csv'),
                'n_asns': 12,
            },
        }

        for test_name, test in tests.items():
            rules = AssociationRegistry(
                global_constraints=test['constraints']
            )
            assert len(rules) >= 3
            for constraint in test['constraints']:
                for rule in rules:
                    yield helpers.check_in_list, \
                        constraint, \
                        rules[rule].GLOBAL_CONSTRAINTS

            pool = AssociationPool.read(test['pool'])
            asns, orphaned = generate(pool, rules)
            yield helpers.check_equal, len(asns), test['n_asns']
