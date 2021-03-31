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

from jwst.lib import engdb_tools
from jwst.lib.tests.engdb_mock import EngDB_Mocker
from jwst.lib import set_telescope_pointing as stp
from jwst import datamodels
from jwst.tests.helpers import word_precision_check

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
    [-1.00444000e-03, 9.99999496e-01, 3.39649146e-06,
     3.38145836e-03, -3.90000000e-14, 9.99994283e-01,
     9.99993778e-01, 1.00444575e-03, -3.38145665e-03]
)
FSMCORR_EXPECTED = np.zeros((2,))
OBSTIME_EXPECTED = STARTTIME

# ########################
# Database access fixtures
# ########################


@pytest.fixture(scope='module')
def method_gscmd_j3pags():
    """Calculate matricies using the GSCMD_J3PAGS method

    This set was derived from the first valid group of enginerring parameters for exposure
    jw00624028002_02101_00001_nrca1 retrieved from the SDP regression tests for Build 7.1.1.
    """
    # setup inputs
    t_pars = stp.TransformParameters()

    t_pars.guide_star_wcs = stp.WCSRef(ra=241.24294932221, dec=70.66165389073196, pa=None)
    t_pars.pointing = stp.Pointing(
        q=np.array([-0.20954692, -0.6177655, -0.44653177, 0.61242575]),
        j2fgs_matrix=np.array([-9.77300013e-04, 3.38988895e-03, 9.99993777e-01,
                               9.99999522e-01, 8.37175385e-09, 9.77305600e-04,
                               3.30458575e-06, 9.99994254e-01, -3.38988734e-03]),
        fsmcorr=np.array([0.00584114, -0.00432878]),
        obstime=Time(1611628160.325, format='unix'),
        gs_commanded=np.array([-22.40031242,  -8.17869377])
    )
    t_pars.siaf = stp.SIAF(v2_ref=120.525464, v3_ref=-527.543132, v3yangle=-0.5627898, vparity=-1,
                           crpix1=1024.5, crpix2=1024.5, cdelt1=0.03113928, cdelt2=0.03132232,
                           vertices_idl=(-32.1682, 32.0906, 31.6586, -31.7234, -32.1683, -32.1904, 32.0823, 31.9456))

    # Calculate the transforms
    transforms = stp.calc_transforms_gscmd_j3pags(t_pars)

    return transforms, t_pars


@pytest.mark.parametrize(
    'matrix, expected',
    [('m_eci2fgs1', np.array([[0.79338298, 0.5312005, 0.29727004],
                              [-0.58752257, 0.7959905, 0.14565831],
                              [-0.15925035, -0.29021568, 0.9436176]])),
     ('m_fgs12sifov', np.array([[9.99761239e-01, -2.18280166e-02, 1.00096493e-03],
                                [2.18291297e-02, 9.99761094e-01, -1.11492579e-03],
                                [-9.76389175e-04, 1.13650978e-03, 9.99998878e-01]])),
     ('m_sifov_fsm_delta', np.array([[1.00000000e+00, 0.00000000e+00, 2.09865089e-08],
                                     [-5.94309657e-16, 1.00000000e+00, 2.83186526e-08],
                                     [-2.09865089e-08, -2.83186526e-08, 1.00000000e+00]])),
     ('m_eci2v', np.array([[-0.16198517, -0.287996, 0.94383214],
                           [0.80585859, 0.51340828, 0.29496417],
                           [-0.56951974, 0.80837506, 0.14891952]])),
     ('m_v2siaf', np.array([[5.59174025e-04, -9.99951603e-01, -9.82234493e-03],
                            [2.56321412e-03, -9.82088099e-03, 9.99948489e-01],
                            [9.99996559e-01, 5.84321994e-04, -2.55759849e-03]])),
     ('m_eci2siaf', np.array([[-0.80031602, -0.5214848, -0.2958849],
                              [-0.57782002, 0.80255299, 0.14843421],
                              [-0.16005712, -0.2897625, 0.94362038]]))
     ]
)
def test_gscmd_j3pags(method_gscmd_j3pags, matrix, expected):
    """Ensure expected calculate of the specified matrix

    Parameters
    ----------
    method_gscmd_j3pags: Transforms, TransformPars
        The calculated transforms and parameters used to do the calculation.

    matrix : str
        The matrix under examination

    expected : numpy.array
        Expected value of the matrix
    """
    transforms, t_pars = method_gscmd_j3pags
    assert np.allclose(getattr(transforms, matrix), expected)


def test_j3pa_at_gs(method_gscmd_j3pags):
    """Ensure J3PA@GS is as expected"""
    transforms, t_pars = method_gscmd_j3pags

    assert np.allclose(t_pars.guide_star_wcs.pa, 297.3522435208429)


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
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '2017-01-01'
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
    model.meta.observation.date = '2017-01-01'

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
    q_exp = np.array([0.62383733, 0.53552715, -0.49252283, 0.28541008])
    j2fgs_exp = np.array([
        -1.00962794e-03, 9.99999464e-01, 3.41404261e-06,
        3.38429719e-03, 2.85793453e-09, 9.99994300e-01,
        9.99993742e-01, 1.00963370e-03, -3.38429548e-03
    ])
    j2fgs_exp = np.array([
        -1.00962794e-03, 3.38429719e-03, 9.99993742e-01,
        9.99999464e-01, 2.85793453e-09, 1.00963370e-03,
        3.41404261e-06, 9.99994300e-01, -3.38429548e-03
    ])
    fsmcorr_exp = np.array([-0.02558673, -0.00200601])
    obstime_exp = Time(1559582740.4880004, format='unix')

    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime,
     gs_commanded) = stp.get_pointing(
         Time('2019-06-03T17:25:40', format='isot').mjd,
         Time('2019-06-03T17:25:56', format='isot').mjd,
    )

    assert np.allclose(q, q_exp)
    assert np.allclose(j2fgs_matrix, j2fgs_exp)
    assert np.allclose(fsmcorr, fsmcorr_exp)
    assert np.isclose(obstime.value, obstime_exp.value)


def test_get_pointing_fail():
    with pytest.raises(Exception):
        q, j2fgs_matrix, fmscorr, obstime, gs_commanded = stp.get_pointing(47892.0, 48256.0)


def test_get_pointing(eng_db_ngas):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime,
     gs_commanded) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert np.isclose(q, Q_EXPECTED).all()
    assert np.isclose(j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
    assert np.isclose(fsmcorr, FSMCORR_EXPECTED).all()
    assert STARTTIME <= obstime <= ENDTIME


def test_logging(eng_db_ngas, caplog):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime,
     gs_commanded) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
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
     obstime,
     gs_commanded) = stp.get_pointing(ZEROTIME_START.mjd, ENDTIME.mjd, reduce_func=stp.first_pointing)
    assert j2fgs_matrix.any()
    (q_desired,
     j2fgs_matrix_desired,
     fsmcorr_desired,
     obstime,
     gs_commanded) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert np.array_equal(q, q_desired)
    assert np.array_equal(j2fgs_matrix, j2fgs_matrix_desired)
    assert np.array_equal(fsmcorr, fsmcorr_desired)


@pytest.mark.skipif(sys.version_info.major < 3,
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


@pytest.mark.skipif(sys.version_info.major < 3,
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


@pytest.mark.skipif(sys.version_info.major < 3,
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


@pytest.mark.skipif(sys.version_info.major < 3,
                    reason="No URI support in sqlite3")
def test_add_wcs_method_original(eng_db_ngas, data_file, siaf_file=siaf_db):
    """Test using the database and the original, pre-JSOCINT-555 algorithms"""
    stp.add_wcs(data_file, siaf_path=siaf_db, method=stp.Methods.ORIGINAL, j2fgs_transpose=False)

    with datamodels.Level1bModel(data_file) as model:
        assert model.meta.visit.pointing_engdb_quality == 'CALCULATED_ORIGINAL'
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


@pytest.mark.skipif(sys.version_info.major < 3,
                    reason="No URI support in sqlite3")
def test_add_wcs_method_gscmd(eng_db_ngas, data_file, siaf_file=siaf_db):
    """Test using the database and the original, pre-JSOCINT-555 algorithms"""
    stp.add_wcs(data_file, siaf_path=siaf_db, method=stp.Methods.GSCMD, j2fgs_transpose=False)

    with datamodels.Level1bModel(data_file) as model:

        assert model.meta.visit.pointing_engdb_quality == 'CALCULATED_GSCMD'
        assert np.isclose(model.meta.pointing.ra_v1, 347.7610269680282)
        assert np.isclose(model.meta.pointing.dec_v1, -86.86190329281591)
        assert np.isclose(model.meta.pointing.pa_v3, 62.07421054339623)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 346.63329946732654)
        assert np.isclose(model.meta.wcsinfo.crval2, -86.95593300483553)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.2628005754440661)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9648501736260754)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9648501736260754)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.2628005754440661)
        assert model.meta.wcsinfo.v2_ref == 200.0
        assert model.meta.wcsinfo.v3_ref == -350.0
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 42.0
        assert np.isclose(model.meta.wcsinfo.ra_ref, 346.63329946732654)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -86.95593300483553)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 63.236303497892706)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 346.238814600 -86.972238983'
                ' 346.401373707 -86.941784876'
                ' 346.968272880 -86.949985434'
                ' 346.816341294 -86.98058802'
            )
        )


@pytest.mark.skipif(sys.version_info.major < 3,
                    reason="No URI support in sqlite3")
def test_add_wcs_method_gscmd_v3pags(eng_db_ngas, data_file, siaf_file=siaf_db):
    """Test using the database and the original, pre-JSOCINT-555 algorithms"""
    stp.add_wcs(data_file, siaf_path=siaf_db, method=stp.Methods.GSCMD_V3PAGS, j2fgs_transpose=False)

    with datamodels.Level1bModel(data_file) as model:

        assert model.meta.visit.pointing_engdb_quality == 'CALCULATED_GSCMD_V3PAGS'
        assert np.isclose(model.meta.pointing.ra_v1, 347.76099011134437)
        assert np.isclose(model.meta.pointing.dec_v1, -86.86190123281543)
        assert np.isclose(model.meta.pointing.pa_v3, 62.073430219473025)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 346.6332873326963)
        assert np.isclose(model.meta.wcsinfo.crval2, -86.95593178595206)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.0555555e-5)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.26278701123754644)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9648538680675108)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9648538680675108)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.26278701123754644)
        assert model.meta.wcsinfo.v2_ref == 200.0
        assert model.meta.wcsinfo.v3_ref == -350.0
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 42.0
        assert np.isclose(model.meta.wcsinfo.ra_ref, 346.6332873326963)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -86.95593178595206)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 63.23549801501118)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 346.238806945 -86.972238059'
                ' 346.401357930 -86.941783830'
                ' 346.968259026 -86.949983965'
                ' 346.816335662 -86.980586665'
            )
        )


@pytest.mark.skipif(sys.version_info.major < 3,
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
