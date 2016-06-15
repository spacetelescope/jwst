from __future__ import absolute_import

from . import helpers

from .. import (AssociationRegistry, AssociationPool, generate)


class TestGenerate():
    pool_file = helpers.t_path('data/jw93060_20150312T160130_pool.csv')

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_generate(self):
        rules = AssociationRegistry()
        pool = AssociationPool.read(self.pool_file)
        (asns, orphaned) = generate(pool, rules)
        assert len(asns) == 14
        assert len(orphaned) == 36
        serialized = asns[0].serialize() # This will test the validity.
