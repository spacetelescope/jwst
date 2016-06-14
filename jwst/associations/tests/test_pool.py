from __future__ import absolute_import

from .helpers import t_path

from .. import AssociationPool


class TestPool():
    pool_file = t_path('data/jw93060_20150312T160130_pool.csv')

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_pool(self):
        pool = AssociationPool.read(self.pool_file)
        assert len(pool) == 636
