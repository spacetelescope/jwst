"""test_associations: Test of general Association functionality."""
from __future__ import absolute_import

from . import helpers

from ..main import Main


class TestMain():

    pools_size = [
        (helpers.t_path('data/jw93060_20150312T160130_pool.csv'), 14),
        (helpers.t_path('data/jw82600_001_20151107T165901_pool.csv'), 11),
    ]

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_script(self):
        for name, number in self.pools_size:
            gs = Main([name, '--dry_run'])
            yield helpers.check_equal, len(gs.associations), number

    def test_asn_candidates(self):
        pool_name = self.pools_size[0][0]
        gs = Main([pool_name, '--dry_run', '-i', '1'])
        assert len(gs.associations) == 6
        gs = Main([pool_name, '--dry_run', '-i', '1', '2'])
        assert len(gs.associations) == 14

    def test_cross_candidate(self):
        pool_name = self.pools_size[1][0]
        gs = Main([pool_name, '--dry_run'])
        assert len(gs.associations) == 11
        gs = Main([pool_name, '--dry_run', '--cross-candidate-only'])
        assert len(gs.associations) == 5
