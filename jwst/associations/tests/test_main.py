"""test_associations: Test of general Association functionality."""
from __future__ import absolute_import
import pytest

from .helpers import full_pool_rules

from ..main import Main


class TestMain(object):

    def test_script(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run'])
        asns = generated.associations
        assert len(asns) == 37
        found_rules = set(
            asn['asn_rule']
            for asn in asns
        )
        assert 'candidate_Asn_Image' in found_rules
        assert 'candidate_Asn_WFSCMB' in found_rules

    def test_asn_candidates(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run', '-i', 'o001'])
        assert len(generated.associations) == 2
        generated = Main([pool_fname, '--dry-run', '-i', 'o001', 'o002'])
        assert len(generated.associations) == 4

    def test_toomanyoptions(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        with pytest.raises(SystemExit):
            generated = Main([
                pool_fname,
                '--dry-run',
                '--discover',
                '--all-candidates',
                '-i', 'o001',
            ])
        with pytest.raises(SystemExit):
            generated = Main([
                pool_fname,
                '--dry-run',
                '--discover',
                '--all-candidates',
            ])
        with pytest.raises(SystemExit):
            generated = Main([
                pool_fname,
                '--dry-run',
                '--discover',
                '-i', 'o001',
            ])
        with pytest.raises(SystemExit):
            generated = Main([
                pool_fname,
                '--dry-run',
                '--all-candidates',
                '-i', 'o001',
            ])

    def test_discovered(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        full = Main([pool_fname, '--dry-run'])
        candidates = Main([pool_fname, '--dry-run', '--all-candidates'])
        discovered = Main([pool_fname, '--dry-run', '--discover'])
        assert len(full.associations) == len(candidates.associations) + len(discovered.associations)
