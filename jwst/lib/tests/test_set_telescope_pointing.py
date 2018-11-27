"""
Test suite for set_telescope_pointing
"""
import copy
import numpy as np
import os
import sys
import pytest
from tempfile import TemporaryDirectory

import requests_mock

from astropy.table import Table
from astropy.time import Time

from .. import engdb_tools
from .engdb_mock import EngDB_Mocker
from .. import set_telescope_pointing as stp
from ... import datamodels
from ...tests.helpers import word_precision_check

# Setup mock engineering service
STARTTIME = Time('2014-01-03')
ENDTIME = Time('2014-01-04')
ZEROTIME_START = Time('2014-01-01')
ZEROTIME_END = Time('2014-01-02')

# Header defaults
TARG_RA = 345.0
TARG_DEC = -87.0
V2_REF = 200.0
V3_REF = -350.0
V3I_YANG = 42.0
VPARITY = -1

# Get the mock DB
db_path = os.path.join(os.path.dirname(__file__), 'data', 'engdb_ngas')
siaf_db = os.path.join(os.path.dirname(__file__), 'data', 'siaf.db')

# Some expected falues
Q_EXPECTED = np.asarray(
    [-0.36915286, 0.33763282, 0.05758533, 0.86395264]
)
J2FGS_MATRIX_EXPECTED = np.asarray(
    [-1.00444000e-03,   3.38145836e-03,   9.99993778e-01,
     9.99999496e-01,  -3.90000000e-14,   1.00444575e-03,
     3.39649146e-06,   9.99994283e-01,  -3.38145665e-03]
)
FSMCORR_EXPECTED = np.zeros((2,))
OBSTIME_EXPECTED = STARTTIME


@pytest.fixture
def eng_db():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_path):
        engdb = engdb_tools.ENGDB_Service()
        yield engdb

@pytest.fixture
def data_file():
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.wcsinfo.v2_ref = V2_REF
    model.meta.wcsinfo.v3_ref = V3_REF
    model.meta.wcsinfo.v3yangle = V3I_YANG
    model.meta.wcsinfo.vparity = VPARITY
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '1/1/2017'

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'fits.fits')
        model.save(file_path)
        yield file_path


def test_get_pointing_fail():
    with pytest.raises(Exception):
        q, j2fgs_matrix, fmscorr, obstime = stp.get_pointing(47892.0, 48256.0)


def test_get_pointing(eng_db):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(STARTTIME, ENDTIME)
    assert np.isclose(q, Q_EXPECTED).all()
    assert np.isclose(j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
    assert np.isclose(fsmcorr, FSMCORR_EXPECTED).all()
    assert STARTTIME <= obstime <= ENDTIME


def test_get_pointing_list(eng_db):
        results = stp.get_pointing(STARTTIME, ENDTIME, result_type='all')
        assert isinstance(results, list)
        assert len(results) > 0
        assert np.isclose(results[0].q, Q_EXPECTED).all()
        assert np.isclose(results[0].j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
        assert np.isclose(results[0].fsmcorr, FSMCORR_EXPECTED).all()
        assert STARTTIME <= results[0].obstime <= ENDTIME


def test_get_pointing_with_zeros(eng_db):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(ZEROTIME_START, ENDTIME)
    assert j2fgs_matrix.any()
    (q_desired,
     j2fgs_matrix_desired,
     fsmcorr_desired,
     obstime) = stp.get_pointing(STARTTIME, ENDTIME)
    assert np.array_equal(q, q_desired)
    assert np.array_equal(j2fgs_matrix, j2fgs_matrix_desired)
    assert np.array_equal(fsmcorr, fsmcorr_desired)


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_default(data_file):
    try:
        stp.add_wcs(data_file, siaf_path=siaf_db)
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    model = datamodels.Level1bModel(data_file)
    assert model.meta.pointing.ra_v1 == TARG_RA
    assert model.meta.pointing.dec_v1 == TARG_DEC
    assert model.meta.pointing.pa_v3 == 0.
    assert model.meta.wcsinfo.crval1 == TARG_RA
    assert model.meta.wcsinfo.crval2 == TARG_DEC
    assert np.isclose(model.meta.wcsinfo.pc1_1, -1.0)
    assert np.isclose(model.meta.wcsinfo.pc1_2, 0.0)
    assert np.isclose(model.meta.wcsinfo.pc2_1, 0.0)
    assert np.isclose(model.meta.wcsinfo.pc2_2, 1.0)
    assert model.meta.wcsinfo.ra_ref == TARG_RA
    assert model.meta.wcsinfo.dec_ref == TARG_DEC
    assert np.isclose(model.meta.wcsinfo.roll_ref, 358.9045979379)
    assert model.meta.wcsinfo.wcsaxes == 2
    assert word_precision_check(
        model.meta.wcsinfo.s_region,
        (
            'POLYGON ICRS'
            ' 345.0516057166881 -86.87312441299257'
            ' 344.61737392066823 -86.85221531104224'
            ' 344.99072891662956 -86.82863042295425'
            ' 345.42745662836063 -86.84915871318734'
        )
    )


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_fsmcorr_v1(data_file):
    """Test with default value using FSM original correction"""
    try:
        stp.add_wcs(data_file, fsmcorr_version='v1', siaf_path=siaf_db)
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    model = datamodels.Level1bModel(data_file)
    assert model.meta.pointing.ra_v1 == TARG_RA
    assert model.meta.pointing.dec_v1 == TARG_DEC
    assert model.meta.pointing.pa_v3 == 0.
    assert model.meta.wcsinfo.crval1 == TARG_RA
    assert model.meta.wcsinfo.crval2 == TARG_DEC
    assert np.isclose(model.meta.wcsinfo.pc1_1, -1.0)
    assert np.isclose(model.meta.wcsinfo.pc1_2, 0.0)
    assert np.isclose(model.meta.wcsinfo.pc2_1, 0.0)
    assert np.isclose(model.meta.wcsinfo.pc2_2, 1.0)
    assert model.meta.wcsinfo.ra_ref == TARG_RA
    assert model.meta.wcsinfo.dec_ref == TARG_DEC
    assert np.isclose(model.meta.wcsinfo.roll_ref, 358.9045979379)
    assert model.meta.wcsinfo.wcsaxes == 2
    assert word_precision_check(
        model.meta.wcsinfo.s_region,
        (
            'POLYGON ICRS'
            ' 345.0516057166881 -86.87312441299257'
            ' 344.61737392066823 -86.85221531104224'
            ' 344.99072891662956 -86.82863042295425'
            ' 345.42745662836063 -86.84915871318734'
        )
    )


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_with_db(eng_db, data_file, siaf_file=siaf_db):
    """Test using the database"""
    stp.add_wcs(data_file, siaf_path=siaf_db)

    model = datamodels.Level1bModel(data_file)
    assert np.isclose(model.meta.pointing.ra_v1, 348.9278669)
    assert np.isclose(model.meta.pointing.dec_v1, -38.749239)
    assert np.isclose(model.meta.pointing.pa_v3, 50.1767077)
    assert np.isclose(model.meta.wcsinfo.crval1, 348.8776709)
    assert np.isclose(model.meta.wcsinfo.crval2, -38.854159)
    assert np.isclose(model.meta.wcsinfo.pc1_1, 0.0385309)
    assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9992574)
    assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9992574)
    assert np.isclose(model.meta.wcsinfo.pc2_2, -0.0385309)
    assert np.isclose(model.meta.wcsinfo.ra_ref, 348.8776709)
    assert np.isclose(model.meta.wcsinfo.dec_ref, -38.854159)
    assert np.isclose(model.meta.wcsinfo.roll_ref, 50.20832726650)
    assert model.meta.wcsinfo.wcsaxes == 2
    assert word_precision_check(
        model.meta.wcsinfo.s_region,
        (
            'POLYGON ICRS'
            ' 349.00694612561705 -38.776964589744054'
            ' 349.0086451128466 -38.74533844552814'
            ' 349.04874980331374 -38.746495669763334'
            ' 349.0474396482846 -38.77812380255898'
        )
    )


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_with_db_fsmcorr_v1(eng_db, data_file):
    """Test using the database with original FSM correction"""
    stp.add_wcs(data_file, fsmcorr_version='v1', siaf_path=siaf_db)

    model = datamodels.Level1bModel(data_file)
    assert np.isclose(model.meta.pointing.ra_v1, 348.9278669)
    assert np.isclose(model.meta.pointing.dec_v1, -38.749239)
    assert np.isclose(model.meta.pointing.pa_v3, 50.1767077)
    assert np.isclose(model.meta.wcsinfo.crval1, 348.8776709)
    assert np.isclose(model.meta.wcsinfo.crval2, -38.854159)
    assert np.isclose(model.meta.wcsinfo.pc1_1, 0.0385309)
    assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9992574)
    assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9992574)
    assert np.isclose(model.meta.wcsinfo.pc2_2, -0.0385309)
    assert np.isclose(model.meta.wcsinfo.ra_ref, 348.8776709)
    assert np.isclose(model.meta.wcsinfo.dec_ref, -38.854159)
    assert np.isclose(model.meta.wcsinfo.roll_ref, 50.20832726650)
    assert model.meta.wcsinfo.wcsaxes == 2
    assert word_precision_check(
        model.meta.wcsinfo.s_region,
        (
            'POLYGON ICRS'
            ' 349.00694612561705 -38.776964589744054'
            ' 349.0086451128466 -38.74533844552814'
            ' 349.04874980331374 -38.746495669763334'
            ' 349.0474396482846 -38.77812380255898'
        )
    )
