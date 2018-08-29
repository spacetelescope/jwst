"""
Test for utc_to_bary
"""
import numpy as np
import os
import pytest
from tempfile import TemporaryDirectory

from .. import utc_to_tdb
from ... import datamodels


@pytest.fixture
def data_file():
    model = datamodels.Level1bModel()
    model.meta.target.ra = 0.
    model.meta.target.dec = 0.
    model.meta.ephemeris.time = 55727.0
    model.meta.ephemeris.spatial_x = -34305.4075983316
    model.meta.ephemeris.spatial_y = 1049528.04998405
    model.meta.ephemeris.spatial_z = 679175.58185602
    model.meta.ephemeris.velocity_x = -0.548663244644384
    model.meta.ephemeris.velocity_y = -0.103904924724239
    model.meta.ephemeris.velocity_z = 0.000982870964178323

    # Assign dummy values to the last three columns.
    model.int_times = \
[(1, 55728.        , 55728.00032119, 55728.00064237, -1., -1., -1.),
 (2, 55728.00064237, 55728.00096356, 55728.00128474, -1., -1., -1.),
 (3, 55728.00128474, 55728.00160593, 55728.00192711, -1., -1., -1.),
 (4, 55728.00192711, 55728.0022483 , 55728.00256948, -1., -1., -1.)]

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'int_times.fits')
        model.save(file_path)
        yield file_path


def test_utc_to_tdb(data_file):
    # This function populates the last three columns in int_times.
    utc_to_tdb.utc_tdb(data_file)

    model = datamodels.open(data_file)
    assert np.isclose(model.meta.ephemeris.time, 55727.0, rtol=1.e-10)
    assert np.isclose(model.meta.ephemeris.spatial_x, -34305.4075983316,
                      rtol=1.e-10)
    assert np.isclose(model.meta.ephemeris.spatial_y, 1049528.04998405,
                      rtol=1.e-10)
    assert np.isclose(model.meta.ephemeris.spatial_z, 679175.58185602,
                      rtol=1.e-10)
    assert np.isclose(model.meta.ephemeris.velocity_x, -0.548663244644384,
                      rtol=1.e-10)
    assert np.isclose(model.meta.ephemeris.velocity_y, -0.103904924724239,
                      rtol=1.e-10)
    assert np.isclose(model.meta.ephemeris.velocity_z, 0.000982870964178323,
                      rtol=1.e-10)
    # These are the last three columns.
    start_tdb = model.int_times['int_start_BJD_TDB']
    mid_tdb = model.int_times['int_mid_BJD_TDB']
    end_tdb = model.int_times['int_end_BJD_TDB']
    # atol = 1.e-6 for values in MJD corresponds to a distance of
    # 1.e-6 * 86400. * 299792.458, which is about 25902 km.
    assert np.allclose(start_tdb, np.array(
        [55728.00016781, 55728.00081024, 55728.00145267, 55728.00209511]),
        rtol=0., atol=1.e-6)
    assert np.allclose(mid_tdb, np.array(
        [55728.00048903, 55728.00113146, 55728.0017739, 55728.00241633]),
        rtol=0., atol=1.e-6)
    assert np.allclose(end_tdb, np.array(
        [55728.00081024, 55728.00145267, 55728.00209511, 55728.00273754]),
        rtol=0., atol=1.e-6)
    model.close()
