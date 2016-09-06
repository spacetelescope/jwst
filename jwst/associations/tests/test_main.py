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
        assert len(asns) == 45
        found_rules = set(
            asn['asn_rule']
            for asn in asns
        )
        assert 'Asn_Image' in found_rules
        assert 'Asn_WFSCMB' in found_rules

    def test_asn_candidates(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run', '-i', 'o001'])
        assert len(generated.associations) == 2
        generated = Main([pool_fname, '--dry-run', '-i', 'o001', 'o002'])
        assert len(generated.associations) == 4

    @pytest.mark.xfail()
    def test_cross_candidate(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run'])
        assert len(generated.associations) == 45
        generated = Main([pool_fname, '--dry-run', '--cross-candidate-only'])
        assert len(generated.associations) == 5
