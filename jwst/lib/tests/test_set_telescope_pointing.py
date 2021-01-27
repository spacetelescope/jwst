"""
Test suite for set_telescope_pointing
"""
import logging
import numpy as np
import os
import sys
import pytest
from tempfile import TemporaryDirectory

from astropy.time import Time

from .. import engdb_tools
from .engdb_mock import EngDB_Mocker
from .. import set_telescope_pointing as stp
from ... import datamodels
from ...tests.helpers import word_precision_check

# Ensure that `set_telescope_pointing` logs.
stp.logger.setLevel(logging.DEBUG)
stp.logger.addHandler(logging.StreamHandler())

# Setup mock engineering service
STARTTIME = Time('2014-01-03')
ENDTIME = Time('2014-01-04')
ZEROTIME_START = Time('2014-01-01')
ZEROTIME_END = Time('2014-01-02')

# Header defaults
TARG_RA = 345.0
TARG_DEC = -87.0

# Get the mock databases
db_ngas_path = os.path.join(os.path.dirname(__file__), 'data', 'engdb_ngas')
db_jw703_path = os.path.join(os.path.dirname(__file__), 'data', 'engdb_jw00703')
siaf_db = os.path.join(os.path.dirname(__file__), 'data', 'siaf.db')

# Some expected falues
Q_EXPECTED = np.asarray(
    [-0.36915286, 0.33763282, 0.05758533, 0.86395264]
)
J2FGS_MATRIX_EXPECTED = np.asarray(
    [-1.00444000e-03, 9.99999496e-01,  3.39649146e-06,
     3.38145836e-03,  -3.90000000e-14, 9.99994283e-01,
     9.99993778e-01,  1.00444575e-03,  -3.38145665e-03]
)
FSMCORR_EXPECTED = np.zeros((2,))
OBSTIME_EXPECTED = STARTTIME

# ########################
# Database access fixtures
# ########################
@pytest.fixture
def eng_db_ngas():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_ngas_path):
        engdb = engdb_tools.ENGDB_Service()
        yield engdb


@pytest.fixture
def eng_db_jw703():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_jw703_path):
        engdb = engdb_tools.ENGDB_Service()
        yield engdb


@pytest.fixture(scope='module')
def data_file():
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '1/1/2017'
    model.meta.exposure.type = "MIR_IMAGE"

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'fits.fits')
        model.save(file_path)
        model.close()
        yield file_path


@pytest.fixture(scope='module')
def data_file_nosiaf():
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.aperture.name = "UNKNOWN"
    model.meta.observation.date = '1/1/2017'

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'fits_nosiaf.fits')
        model.save(file_path)
        model.close()
        yield file_path


def test_change_engdb_url():
    """Test changing the engineering database by call for success.

    The given time and database should not find any values.
    """
    with pytest.raises(ValueError):
        stp.get_pointing(
            STARTTIME.mjd,
            ENDTIME.mjd,
            engdb_url=engdb_tools.ENGDB_BASE_URL
        )


def test_change_engdb_url_fail():
    """Test changing the engineering database by call"""
    with pytest.raises(Exception):
        stp.get_pointing(
            Time('2019-06-03T17:25:40', format='isot').mjd,
            Time('2019-06-03T17:25:56', format='isot').mjd,
            engdb_url='http://nonexistant.fake'
        )


def test_strict_pointing(data_file, eng_db_jw703):
    """Test failure on strict pointing"""
    with pytest.raises(ValueError):
        stp.add_wcs(data_file, siaf_path=siaf_db, tolerance=0)


def test_pointing_averaging(eng_db_jw703):
    """Ensure that the averaging works."""
    q_exp = np.array([ 0.62383733,  0.53552715, -0.49252283,  0.28541008])
    j2fgs_exp = np.array([
        -1.00962794e-03,  9.99999464e-01,  3.41404261e-06,
        3.38429719e-03, 2.85793453e-09,  9.99994300e-01,
        9.99993742e-01,  1.00963370e-03, -3.38429548e-03
    ])
    j2fgs_exp = np.array([
        -1.00962794e-03, 3.38429719e-03, 9.99993742e-01,
        9.99999464e-01,  2.85793453e-09, 1.00963370e-03,
        3.41404261e-06,  9.99994300e-01, -3.38429548e-03
    ])
    fsmcorr_exp = np.array([-0.02558673, -0.00200601])
    obstime_exp = Time(1559582740.4880004, format='unix')

    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(
         Time('2019-06-03T17:25:40', format='isot').mjd,
         Time('2019-06-03T17:25:56', format='isot').mjd,
     )

    assert np.allclose(q, q_exp)
    assert np.allclose(j2fgs_matrix, j2fgs_exp)
    assert np.allclose(fsmcorr, fsmcorr_exp)
    assert np.isclose(obstime.value, obstime_exp.value)


def test_get_pointing_fail():
    with pytest.raises(Exception):
        q, j2fgs_matrix, fmscorr, obstime = stp.get_pointing(47892.0, 48256.0)


def test_get_pointing(eng_db_ngas):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert np.isclose(q, Q_EXPECTED).all()
    assert np.isclose(j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
    assert np.isclose(fsmcorr, FSMCORR_EXPECTED).all()
    assert STARTTIME <= obstime <= ENDTIME


def test_logging(eng_db_ngas, caplog):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert 'Determining pointing between observations times' in caplog.text
    assert 'Telemetry search tolerance' in caplog.text
    assert 'Reduction function' in caplog.text
    assert 'Querying engineering DB' in caplog.text


def test_get_pointing_list(eng_db_ngas):
    results = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, reduce_func=stp.all_pointings)
    assert isinstance(results, list)
    assert len(results) > 0
    assert np.isclose(results[0].q, Q_EXPECTED).all()
    assert np.isclose(results[0].j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
    assert np.isclose(results[0].fsmcorr, FSMCORR_EXPECTED).all()
    assert STARTTIME <= results[0].obstime <= ENDTIME


def test_get_pointing_with_zeros(eng_db_ngas):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(ZEROTIME_START.mjd, ENDTIME.mjd, reduce_func=stp.first_pointing)
    assert j2fgs_matrix.any()
    (q_desired,
     j2fgs_matrix_desired,
     fsmcorr_desired,
     obstime) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert np.array_equal(q, q_desired)
    assert np.array_equal(j2fgs_matrix, j2fgs_matrix_desired)
    assert np.array_equal(fsmcorr, fsmcorr_desired)


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_default(data_file):
    """Handle when no pointing exists and the default is used."""
    try:
        stp.add_wcs(
            data_file, siaf_path=siaf_db, tolerance=0, allow_default=True
        )
    except ValueError:
        pass  # This is what we want for the test.
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    with datamodels.Level1bModel(data_file) as model:
        assert model.meta.pointing.ra_v1 == TARG_RA
        assert model.meta.pointing.dec_v1 == TARG_DEC
        assert model.meta.pointing.pa_v3 == 0.
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert model.meta.wcsinfo.crval1 == TARG_RA
        assert model.meta.wcsinfo.crval2 == TARG_DEC
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.pc1_1, -0.7558009243361943)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.654801468211972)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.654801468211972)
        assert np.isclose(model.meta.wcsinfo.pc2_2, 0.7558009243361943)
        assert model.meta.wcsinfo.v2_ref == 200.0
        assert model.meta.wcsinfo.v3_ref == -350.0
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 42.0
        assert model.meta.wcsinfo.ra_ref == TARG_RA
        assert model.meta.wcsinfo.dec_ref == TARG_DEC
        assert np.isclose(model.meta.wcsinfo.roll_ref, 358.9045979379)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 345.11054995209815 -87.02586884935684'
                ' 344.6537904121288 -87.00498014679253'
                ' 345.04569816117015 -86.98138111042982'
                ' 345.50498899320183 -87.00187988107017'
            )
        )


def test_add_wcs_default_nosiaf(data_file_nosiaf, caplog):
    """Handle when no pointing exists and the default is used and no SIAF specified."""
    with pytest.raises(ValueError):
        stp.add_wcs(
            data_file_nosiaf, siaf_path=siaf_db, tolerance=0, allow_default=True
        )


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_fsmcorr_v1(data_file):
    """Test with default value using FSM original correction"""
    try:
        stp.add_wcs(
            data_file, fsmcorr_version='v1', siaf_path=siaf_db, tolerance=0, allow_default=True
        )
    except ValueError:
        pass  # This is what we want for the test.
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    with  datamodels.Level1bModel(data_file) as model:
        assert model.meta.pointing.ra_v1 == TARG_RA
        assert model.meta.pointing.dec_v1 == TARG_DEC
        assert model.meta.pointing.pa_v3 == 0.
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert model.meta.wcsinfo.crval1 == TARG_RA
        assert model.meta.wcsinfo.crval2 == TARG_DEC
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.pc1_1, -0.7558009243361943)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.654801468211972)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.654801468211972)
        assert np.isclose(model.meta.wcsinfo.pc2_2, 0.7558009243361943)
        assert model.meta.wcsinfo.v2_ref == 200.0
        assert model.meta.wcsinfo.v3_ref == -350.0
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 42.0
        assert model.meta.wcsinfo.ra_ref == TARG_RA
        assert model.meta.wcsinfo.dec_ref == TARG_DEC
        assert np.isclose(model.meta.wcsinfo.roll_ref, 358.9045979379)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 345.11054995209815 -87.02586884935684'
                ' 344.6537904121288 -87.00498014679253'
                ' 345.04569816117015 -86.98138111042982'
                ' 345.50498899320183 -87.00187988107017'
            )
        )


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_with_db(eng_db_ngas, data_file, siaf_file=siaf_db):
    """Test using the database"""
    stp.add_wcs(data_file, siaf_path=siaf_db, j2fgs_transpose=False)

    with datamodels.Level1bModel(data_file) as model:
        assert np.isclose(model.meta.pointing.ra_v1, 348.9278669)
        assert np.isclose(model.meta.pointing.dec_v1, -38.749239)
        assert np.isclose(model.meta.pointing.pa_v3, 50.1767077)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 348.8776709)
        assert np.isclose(model.meta.wcsinfo.crval2, -38.854159)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.03853303979862607)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9992573266400789)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9992573266400789)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.03853303979862607)
        assert model.meta.wcsinfo.v2_ref == 200.0
        assert model.meta.wcsinfo.v3_ref == -350.0
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 42.0
        assert np.isclose(model.meta.wcsinfo.ra_ref, 348.8776709)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -38.854159)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 50.20832726650)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 348.8563379013152 -38.874810886750495'
                ' 348.85810582665334 -38.84318773861823'
                ' 348.8982592685148 -38.84439628911871'
                ' 348.89688051688233 -38.876020020321164'
            )
        )


@pytest.mark.skipif(sys.version_info.major<3,
                    reason="No URI support in sqlite3")
def test_add_wcs_with_db_fsmcorr_v1(eng_db_ngas, data_file):
    """Test using the database with original FSM correction"""
    stp.add_wcs(data_file, fsmcorr_version='v1', siaf_path=siaf_db, j2fgs_transpose=False)

    with datamodels.Level1bModel(data_file) as model:
        assert np.isclose(model.meta.pointing.ra_v1, 348.9278669)
        assert np.isclose(model.meta.pointing.dec_v1, -38.749239)
        assert np.isclose(model.meta.pointing.pa_v3, 50.1767077)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 348.8776709)
        assert np.isclose(model.meta.wcsinfo.crval2, -38.854159)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.03853303979862607)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9992573266400789)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9992573266400789)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.03853303979862607)
        assert model.meta.wcsinfo.v2_ref == 200.0
        assert model.meta.wcsinfo.v3_ref == -350.0
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 42.0
        assert np.isclose(model.meta.wcsinfo.ra_ref, 348.8776709)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -38.854159)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 50.20832726650)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 348.8563379013152 -38.874810886750495'
                ' 348.85810582665334 -38.84318773861823'
                ' 348.8982592685148 -38.84439628911871'
                ' 348.89688051688233 -38.876020020321164'
            )
        )


def test_default_siaf_values(eng_db_ngas, data_file_nosiaf):
    """
    Test that FITS WCS default values were set.
    """
    with datamodels.Level1bModel(data_file_nosiaf) as model:
        model.meta.exposure.start_time = STARTTIME.mjd
        model.meta.exposure.end_time = ENDTIME.mjd
        model.meta.target.ra = TARG_RA
        model.meta.target.dec = TARG_DEC
        model.meta.aperture.name = "MIRIM_TAFULL"
        model.meta.observation.date = '1/1/2017'
        model.meta.exposure.type = "MIR_IMAGE"
        stp.update_wcs(model, siaf_path=siaf_db, allow_default=False)
        assert model.meta.wcsinfo.crpix1 == 0
        assert model.meta.wcsinfo.crpix2 == 0
        assert model.meta.wcsinfo.cdelt1 == 1
        assert model.meta.wcsinfo.cdelt2 == 1


def test_tsgrism_siaf_values(eng_db_ngas, data_file_nosiaf):
    """
    Test that FITS WCS default values were set.
    """
    with datamodels.Level1bModel(data_file_nosiaf) as model:
        model.meta.exposure.start_time = STARTTIME.mjd
        model.meta.exposure.end_time = ENDTIME.mjd
        model.meta.aperture.name = "NRCA5_GRISM256_F444W"
        model.meta.observation.date = '1/1/2017'
        model.meta.exposure.type = "NRC_TSGRISM"
        model.meta.visit.tsovisit = True
        stp.update_wcs(model, siaf_path=siaf_db)
        assert model.meta.wcsinfo.siaf_xref_sci == 887
        assert model.meta.wcsinfo.siaf_yref_sci == 35
