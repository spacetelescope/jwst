from nose.tools import assert_raises

from jwst_tools.associations.pool import AssociationPool


class TestPool():
    pool_file = 'tests/data/jw93060_20150312T160130_pool.csv'

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_pool(self):
        pool = AssociationPool.read(self.pool_file)
        assert len(pool) == 636
