import os
from jwst.associations.tests.helpers import t_path
from jwst.associations import AssociationPool

POOL_FILE = t_path('data/jw93060_20150312T160130_pool.csv')


def test_pool(tmp_path_factory):
    pool = AssociationPool.read(POOL_FILE)
    assert len(pool) == 636

    tmp_path = tmp_path_factory.mktemp(__name__)
    tmp_pool = tmp_path / 'tmp_pool.csv'
    pool.write(tmp_pool)

    roundtrip = AssociationPool.read(tmp_pool)
    assert len(pool) == len(roundtrip)
    assert set(pool.colnames) == set(roundtrip.colnames)
