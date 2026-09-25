from astropy.utils.data import get_pkg_data_filename

from jwst.associations import AssociationPool


def test_pool(tmp_path):
    pool = AssociationPool.read(
        get_pkg_data_filename(
            "data/jw93060_20150312T160130_pool.csv", package="jwst.associations.tests"
        )
    )
    assert len(pool) == 636

    tmp_pool = tmp_path / "tmp_pool.csv"
    pool.write(tmp_pool)

    roundtrip = AssociationPool.read(tmp_pool)
    assert len(pool) == len(roundtrip)
    assert set(pool.colnames) == set(roundtrip.colnames)
