"""
Test suite for set_telescope_pointing
"""
import logging
import numpy as np
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import pytest

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
DATA_PATH = Path(__file__).parent / 'data'
db_ngas_path = DATA_PATH  / 'engdb_ngas'
db_jw703_path = DATA_PATH / 'engdb_jw00703'
siaf_db = DATA_PATH / 'siaf.db'

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


# Expected transforms
TFORMS_FULL = stp.Transforms.from_asdf(DATA_PATH / 'tforms_full.asdf')
TFORMS_FULLVA = stp.Transforms.from_asdf(DATA_PATH / 'tforms_fullva.asdf')
TFORMS_GSCMD_J3PAGS = stp.Transforms.from_asdf(DATA_PATH / 'tforms_gscmd_j3pags.asdf')


@pytest.fixture(scope='module',
                params=[method for method in stp.Methods])
def calc_transforms(request, tmp_path_factory):
    """Calculate matrices for specified method method

    This set was derived from the first valid group of engineering parameters for exposure
    jw00624028002_02101_00001_nrca1 retrieved from the SDP regression tests for Build 7.7.1.
    """
    # setup inputs
    t_pars = stp.TransformParameters()
    t_pars.method = request.param

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
    t_pars.siaf_path = siaf_db

    # Calculate the transforms
    transforms = stp.calc_transforms(t_pars)

    # Save transforms for later examination
    transforms.write_to_asdf(tmp_path_factory.mktemp('transforms') / f'tforms_{request.param}.asdf')

    return transforms, t_pars


@pytest.mark.parametrize(
    'method',
    [method for method in stp.Methods]
)
def test_method_string(method):
    """Ensure that the value of the method is the string representation"""
    assert f'{method}' == method.value


def test_override_calc_wcs(method_fullva):
    """Test matrix override in the full calculation"""
    transforms, t_pars = method_fullva
    wcsinfo, vinfo, _ = stp.calc_wcs(t_pars)

    override = stp.Transforms(m_eci2j=np.array([[0.80583682, 0.51339893, 0.29503999],
                                                [-0.56953229, 0.8083677, 0.14891175],
                                                [-0.16204967, -0.28803339, 0.94380971]]))
    t_pars.override_transforms = override
    wcsinfo_new, vinfo_new, transforms_new = stp.calc_wcs(t_pars)

    assert vinfo_new != vinfo
    assert vinfo_new == stp.WCSRef(ra=32.5101542644757, dec=17.16214932779597, pa=352.2869379284309)


@pytest.mark.parametrize(
    'attribute, expected',
    [('m_eci2j', 'overridden'), ('m_j2fgs1', 'untouched')]
)
def test_override(attribute, expected):
    """Test overriding of Transforms attributes"""
    overrides = stp.Transforms(m_eci2j='overridden')
    to_override = stp.Transforms(m_eci2j='original', m_j2fgs1='untouched', override=overrides)

    assert getattr(to_override, attribute) == expected


def test_transform_serialize(calc_transforms, tmp_path):
    """Test serialization of Transforms"""
    transforms, t_pars = calc_transforms

    path = tmp_path / 'transforms.asdf'
    transforms.write_to_asdf(path)
    from_asdf = stp.Transforms.from_asdf(path)

    assert isinstance(from_asdf, stp.Transforms)
    assert str(transforms) == str(from_asdf)


def test_methods(fixture, matrix, expected, request):
    """Ensure expected calculate of the specified matrix

    Parameters
    ----------
    fixture: str
        The fixture to retrieve calculated transforms and parameters used to do the calculation.

    matrix : str
        The matrix under examination

    expected : numpy.array
        Expected value of the matrix
    """
    transforms, t_pars = request.getfixturevalue(fixture)
    value = getattr(transforms, matrix)
    if expected is None:
        assert value is None
    else:
        assert np.allclose(value, expected)


def test_j3pa_at_gs(method_gscmd_j3pags):
    """Ensure J3PA@GS is as expected"""
    transforms, t_pars = method_gscmd_j3pags

    assert np.allclose(t_pars.guide_star_wcs.pa, 297.3522435208429)


@pytest.fixture(scope='module')
def eng_db_ngas():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_ngas_path):
        engdb = engdb_tools.ENGDB_Service()
        yield engdb


@pytest.fixture(scope='module')
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
    model.meta.ephemeris.velocity_x = -25.021
    model.meta.ephemeris.velocity_y = -16.507
    model.meta.ephemeris.velocity_z = -7.187

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'fits.fits')
        model.save(file_path)
        model.close()
        yield file_path


@pytest.fixture
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
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.067124166666667e-05)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.090061944444444e-05)
        assert np.isclose(model.meta.wcsinfo.pc1_1, -0.9918437873477998)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.12745941118478668)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.12745941118478668)
        assert np.isclose(model.meta.wcsinfo.pc2_2, 0.9918437873477998)
        assert model.meta.wcsinfo.v2_ref == -453.559116
        assert model.meta.wcsinfo.v3_ref == -373.814447
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 4.83425324
        assert model.meta.wcsinfo.ra_ref == TARG_RA
        assert model.meta.wcsinfo.dec_ref == TARG_DEC
        assert np.isclose(model.meta.wcsinfo.roll_ref, 2.4885527140636143)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 345.367907569 -87.018095946'
                ' 344.762655544 -87.014181559'
                ' 344.844131586 -86.983147828'
                ' 345.436760816 -86.987021446'
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
def test_add_wcs_with_db(eng_db_ngas, data_file):
    """Test using the database"""
    stp.add_wcs(data_file, siaf_path=siaf_db)

    with datamodels.Level1bModel(data_file) as model:
        assert np.isclose(model.meta.pointing.ra_v1, 115.9149987243166)
        assert np.isclose(model.meta.pointing.dec_v1, -37.10479476639234)
        assert np.isclose(model.meta.pointing.pa_v3, 218.59009356876257)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 116.11967778822046)
        assert np.isclose(model.meta.wcsinfo.crval2, -37.10204122234404)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.067124166666667e-05)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.090061944444444e-05)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.7277635275872516)
        assert np.isclose(model.meta.wcsinfo.pc1_2, -0.685828147507639)
        assert np.isclose(model.meta.wcsinfo.pc2_1, -0.685828147507639)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.7277635275872516)
        assert model.meta.wcsinfo.v2_ref == -453.559116
        assert model.meta.wcsinfo.v3_ref == -373.814447
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 4.83425324
        assert np.isclose(model.meta.wcsinfo.ra_ref, 116.11967778822046)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -37.10204122234404)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 218.4665178683768)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 116.113635378 -37.076103949'
                ' 116.142685295 -37.097804891'
                ' 116.115595346 -37.120488860'
                ' 116.086890558 -37.098982235'
            )
        )


@pytest.mark.skipif(sys.version_info.major < 3,
                    reason="No URI support in sqlite3")
def test_add_wcs_method_gscmd(eng_db_ngas, data_file):
    """Test using the database and the original, pre-JSOCINT-555 algorithms"""

    # Calculate
    stp.add_wcs(data_file, siaf_path=siaf_db, method=stp.Methods.GSCMD)

    # Test
    with datamodels.Level1bModel(data_file) as model:

        assert model.meta.visit.pointing_engdb_quality == 'CALCULATED_GSCMD'
        assert np.isclose(model.meta.pointing.ra_v1, 347.7610269680282)
        assert np.isclose(model.meta.pointing.dec_v1, -86.86190329281591)
        assert np.isclose(model.meta.pointing.pa_v3, 62.07421054339623)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 345.0631069614913)
        assert np.isclose(model.meta.wcsinfo.crval2, -86.79567097531229)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.067124166666667e-05)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.090061944444444e-05)
        assert np.isclose(model.meta.wcsinfo.pc1_1, -0.3494308313783168)
        assert np.isclose(model.meta.wcsinfo.pc1_2, 0.9369621625670155)
        assert np.isclose(model.meta.wcsinfo.pc2_1, 0.9369621625670155)
        assert np.isclose(model.meta.wcsinfo.pc2_2, 0.3494308313783168)
        assert model.meta.wcsinfo.v2_ref == -453.559116
        assert model.meta.wcsinfo.v3_ref == -373.814447
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 4.83425324
        assert np.isclose(model.meta.wcsinfo.ra_ref, 345.0631069614913)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -86.79567097531229)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 64.71324058202622)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 344.934234572 -86.821057863'
                ' 344.735601211 -86.791300754'
                ' 345.260360839 -86.780545799'
                ' 345.460312293 -86.809898942'
            )
        )


def test_add_wcs_method_full_nosiafdb(eng_db_ngas, data_file):
    """Test using the database and the original, post-JSOCINT-555 quaternion-based algorithm"""
    # Only run if `pysiaf` is installed.
    pytest.importorskip('pysiaf')

    # Calculate
    stp.add_wcs(data_file, method=stp.Methods.FULL)

    # Test
    with datamodels.Level1bModel(data_file) as model:

        print(f'*****\n{model.meta.instance}\n*****')

        assert model.meta.visit.pointing_engdb_quality == 'CALCULATED_FULL'
        assert np.isclose(model.meta.pointing.ra_v1, 115.9149987243166)
        assert np.isclose(model.meta.pointing.dec_v1, -37.10493815094143)
        assert np.isclose(model.meta.pointing.pa_v3, 218.59229460472764)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 116.11989889855879)
        assert np.isclose(model.meta.wcsinfo.crval2, -37.10219087673784)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.067124166666667e-05)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.090061944444444e-05)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.7277371898873619)
        assert np.isclose(model.meta.wcsinfo.pc1_2, -0.6858560945671079)
        assert np.isclose(model.meta.wcsinfo.pc2_1, -0.6858560945671079)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.7277371898873619)
        assert model.meta.wcsinfo.v2_ref == -453.559116
        assert model.meta.wcsinfo.v3_ref == -373.814447
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 4.83425324
        assert np.isclose(model.meta.wcsinfo.ra_ref, 116.11989889855879)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -37.10219087673784)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 218.46871814012928)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 116.113857725 -37.076253418'
                ' 116.142906655 -37.097955250'
                ' 116.115815560 -37.120638389'
                ' 116.087111751 -37.099130885'
            )
        )


@pytest.mark.skipif(sys.version_info.major < 3,
                    reason="No URI support in sqlite3")
def test_add_wcs_method_full_siafdb(eng_db_ngas, data_file):
    """Test using the database and the original, post-JSOCINT-555 quaternion-based algorithm"""

    # Calculate
    stp.add_wcs(data_file, siaf_path=siaf_db, method=stp.Methods.FULL)

    # Test
    with datamodels.Level1bModel(data_file) as model:

        print(f'*****\n{model.meta.instance}\n*****')

        assert model.meta.visit.pointing_engdb_quality == 'CALCULATED_FULL'
        assert np.isclose(model.meta.pointing.ra_v1, 115.91521930595907)
        assert np.isclose(model.meta.pointing.dec_v1, -37.10493815094143)
        assert np.isclose(model.meta.pointing.pa_v3, 218.59229460472764)
        assert model.meta.wcsinfo.wcsaxes == 2
        assert model.meta.wcsinfo.crpix1 == 693.5
        assert model.meta.wcsinfo.crpix2 == 512.5
        assert np.isclose(model.meta.wcsinfo.crval1, 116.11989889855879)
        assert np.isclose(model.meta.wcsinfo.crval2, -37.10219087673784)
        assert model.meta.wcsinfo.ctype1 == "RA---TAN"
        assert model.meta.wcsinfo.ctype2 == "DEC--TAN"
        assert model.meta.wcsinfo.cunit1 == 'deg'
        assert model.meta.wcsinfo.cunit2 == 'deg'
        assert np.isclose(model.meta.wcsinfo.cdelt1, 3.067124166666667e-05)
        assert np.isclose(model.meta.wcsinfo.cdelt2, 3.090061944444444e-05)
        assert np.isclose(model.meta.wcsinfo.pc1_1, 0.7277371898873619)
        assert np.isclose(model.meta.wcsinfo.pc1_2, -0.6858560945671079)
        assert np.isclose(model.meta.wcsinfo.pc2_1, -0.6858560945671079)
        assert np.isclose(model.meta.wcsinfo.pc2_2, -0.7277371898873619)
        assert model.meta.wcsinfo.v2_ref == -453.559116
        assert model.meta.wcsinfo.v3_ref == -373.814447
        assert model.meta.wcsinfo.vparity == -1
        assert model.meta.wcsinfo.v3yangle == 4.83425324
        assert np.isclose(model.meta.wcsinfo.ra_ref, 116.11989889855879)
        assert np.isclose(model.meta.wcsinfo.dec_ref, -37.10219087673784)
        assert np.isclose(model.meta.wcsinfo.roll_ref, 218.46871814012928)
        assert word_precision_check(
            model.meta.wcsinfo.s_region,
            (
                'POLYGON ICRS'
                ' 116.113857725 -37.076253418'
                ' 116.142906655 -37.097955250'
                ' 116.115815560 -37.120638389'
                ' 116.087111751 -37.099130885'
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
        model.meta.observation.date = '2017-01-01'
        model.meta.exposure.type = "MIR_IMAGE"
        stp.update_wcs(model, siaf_path=siaf_db, allow_default=False)
        assert model.meta.wcsinfo.crpix1 == 24.5
        assert model.meta.wcsinfo.crpix2 == 24.5
        assert model.meta.wcsinfo.cdelt1 == 3.067124166666667e-05
        assert model.meta.wcsinfo.cdelt2 == 3.090061944444444e-05


def test_tsgrism_siaf_values(eng_db_ngas, data_file_nosiaf):
    """
    Test that FITS WCS default values were set.
    """
    with datamodels.Level1bModel(data_file_nosiaf) as model:
        model.meta.exposure.start_time = STARTTIME.mjd
        model.meta.exposure.end_time = ENDTIME.mjd
        model.meta.aperture.name = "NRCA5_GRISM256_F444W"
        model.meta.observation.date = '2017-01-01'
        model.meta.exposure.type = "NRC_TSGRISM"
        model.meta.visit.tsovisit = True
        stp.update_wcs(model, siaf_path=siaf_db)
        assert model.meta.wcsinfo.siaf_xref_sci == 952
        assert model.meta.wcsinfo.siaf_yref_sci == 35
