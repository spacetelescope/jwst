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
        assert len(asns) == 44
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

    def test_all_candidates(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run', '--all-candidates'])
        assert len(generated.associations) == 2

    @pytest.mark.xfail()
    def test_discover(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run', '--discover'])
        assert len(generated.associations) == 2


    @pytest.mark.xfail()
    def test_cross_candidate(self, full_pool_rules):
        pool, rules, pool_fname = full_pool_rules

        generated = Main([pool_fname, '--dry-run'])
        assert len(generated.associations) == 44
        generated = Main([pool_fname, '--dry-run', '--cross-candidate-only'])
        assert len(generated.associations) == 5
