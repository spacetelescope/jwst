from jwst.associations.association import AssociationRegistry
from jwst.associations.pool import AssociationPool
from jwst.associations.generate import generate


class TestGenerate():
    pool_file = 'tests/data/jw93060_20150312T160130_pool.csv'

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
        json = asns[0].to_json() # This will test the validity.
