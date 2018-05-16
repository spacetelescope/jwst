from .helpers import t_path

from .. import AssociationPool

POOL_FILE = t_path('data/jw93060_20150312T160130_pool.csv')


def test_pool(tmpdir):
    pool = AssociationPool.read(POOL_FILE)
    assert len(pool) == 636

    tmp_pool = str(tmpdir.mkdir(__name__).join('tmp_pool.csv'))
    pool.write(tmp_pool)

    roundtrip = AssociationPool.read(tmp_pool)
    assert len(pool) == len(roundtrip)
    assert set(pool.colnames) == set(roundtrip.colnames)
