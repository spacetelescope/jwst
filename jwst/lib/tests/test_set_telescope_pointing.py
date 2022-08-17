"""
Test suite for set_telescope_pointing
"""
import logging
import numpy as np
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

import pytest

from astropy.io import fits
from astropy.time import Time

from jwst import datamodels
from jwst.lib import engdb_mast
from jwst.lib import engdb_tools
from jwst.lib import set_telescope_pointing as stp
from jwst.lib import siafdb
from jwst.lib.tests.engdb_mock import EngDB_Mocker
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
db_ngas_path = DATA_PATH / 'engdb_ngas'
db_jw703_path = DATA_PATH / 'engdb_jw00703'
siaf_path = DATA_PATH / 'siaf.db'

# Some expected values
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

# Meta attributes for test comparisons
METAS_EQUALITY = ['meta.visit.engdb_pointing_quality',
                  'meta.wcsinfo.wcsaxes',
                  'meta.wcsinfo.ctype1',
                  'meta.wcsinfo.ctype2',
                  'meta.wcsinfo.cunit1',
                  'meta.wcsinfo.cunit2',
                  'meta.wcsinfo.vparity',
                  ]
METAS_ISCLOSE = ['meta.wcsinfo.crpix1',
                 'meta.wcsinfo.crpix2',
                 'meta.wcsinfo.crval1',
                 'meta.wcsinfo.crval2',
                 'meta.wcsinfo.cdelt1',
                 'meta.wcsinfo.cdelt2',
                 'meta.wcsinfo.pc1_1',
                 'meta.wcsinfo.pc1_2',
                 'meta.wcsinfo.pc2_1',
                 'meta.wcsinfo.pc2_2',
                 'meta.wcsinfo.roll_ref',
                 'meta.wcsinfo.v2_ref',
                 'meta.wcsinfo.v3_ref',
                 'meta.wcsinfo.v3yangle',
                 'meta.wcsinfo.ra_ref',
                 'meta.wcsinfo.dec_ref',
                 'meta.pointing.ra_v1',
                 'meta.pointing.dec_v1',
                 'meta.pointing.pa_v3',
                 ]


@pytest.fixture(params=[('good_model', True), ('bad_model', False), ('fits_nomodel', False)])
def file_case(request, tmp_path):
    """Generate files with different model states"""
    case, allow = request.param

    if case == 'good_model':
        # Make a model that will always succeed
        model = datamodels.Level1bModel((10, 10, 10, 10))
        path = tmp_path / 'level1bmodel.fits'
        model.save(path)
    elif case == 'bad_model':
        # Make a model that will fail if not allowed
        model = datamodels.IFUCubeModel((10, 10, 10))
        path = tmp_path / 'image.fits'
        model.save(path)
    elif case == 'fits_nomodel':
        # Create just a plain anything FITS
        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        path = tmp_path / 'empty.fits'
        hdul.writeto(path)
    else:
        assert False, f'Cannot produce a file for {case}'

    return path, allow


@pytest.mark.parametrize('allow_any_file', [True, False])
def test_allow_any_file(file_case, allow_any_file):
    """Test various files against whether they should be allowed or not

    Parameters
    ----------
    file_case : (Path-like, allow)
        File to test and whether it should always be allowed.
        If not `allow`, the file should be usable only when `allow_any_file`.

    allow_any_file : bool
        Value of `allow_any_file` to try
    """
    path, allow = file_case
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        if not allow and not allow_any_file:
            with pytest.raises(TypeError):
                stp.add_wcs(path, siaf_path=siaf_path, allow_any_file=allow_any_file, dry_run=True)
        else:
            # Expected error when trying to actually add the wcs.
            # The provided files do not have sufficient info to do the calculations.
            with pytest.raises(AttributeError):
                stp.add_wcs(path, siaf_path=siaf_path, allow_any_file=allow_any_file, dry_run=True)


@pytest.mark.parametrize(
    'method',
    [method for method in stp.Methods]
)
def test_method_string(method):
    """Ensure that the value of the method is the string representation"""
    assert f'{method}' == method.value


def test_override_calc_wcs():
    """Test matrix override in the full calculation"""
    t_pars = make_t_pars()
    t_pars.method = stp.Methods.OPS_TR_202111
    wcsinfo, vinfo, _ = stp.calc_wcs(t_pars)

    override = stp.Transforms(m_eci2j=np.array([[0.80583682, 0.51339893, 0.29503999],
                                                [-0.56953229, 0.8083677, 0.14891175],
                                                [-0.16204967, -0.28803339, 0.94380971]]))
    t_pars.override_transforms = override
    wcsinfo_new, vinfo_new, transforms_new = stp.calc_wcs(t_pars)

    assert vinfo_new != vinfo
    assert all(np.isclose(vinfo_new,
                          stp.WCSRef(ra=32.50407337171124, dec=17.161233048951043, pa=352.28553015159287)))


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


@pytest.mark.parametrize('matrix', [matrix for matrix in stp.Transforms()._fields])
def test_methods(calc_transforms, matrix):
    """Ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare
    """
    _test_methods(calc_transforms, matrix)


def test_coarse_202111_fgsid(calc_coarse_202111_fgsid):
    """Test COARSE_202111 to ensure correct FGS id is used.

    If an FGS is specifically used for science, the other FGS is the guider.
    For all other instruments, the assumption is that FGS1 is the guider, but
    that can be overridden.

    This tests that the correct FGS id was used in the calculations
    """
    transforms, t_pars, truth_ext, fgs_expected = calc_coarse_202111_fgsid

    fgsid_used = t_pars.fgsid
    assert fgsid_used == fgs_expected


@pytest.mark.parametrize('matrix', [matrix for matrix in stp.Transforms()._fields])
def test_coarse_202111_fgsid_matrices(calc_coarse_202111_fgsid, matrix):
    """Test COARSE_202111 matrices using various FGS settings

    If an FGS is specifically used for science, the other FGS is the guider.
    For all other instruments, the assumption is that FGS1 is the guider, but
    that can be overridden.

    This tests that the appropriate transformations were calculated.
    """
    transforms, t_pars, truth_ext, fgs_expected = calc_coarse_202111_fgsid
    _test_methods((transforms, t_pars), matrix, truth_ext=truth_ext)


def test_change_engdb_url():
    """Test changing the engineering database by call for success.

    The given time and database should not find any values.
    """
    with pytest.raises(ValueError):
        stp.get_pointing(
            STARTTIME.mjd,
            ENDTIME.mjd,
            engdb_url=engdb_mast.MAST_BASE_URL
        )


def test_change_engdb_url_fail():
    """Test changing the engineering database by call"""
    with pytest.raises(Exception):
        stp.get_pointing(
            Time('2019-06-03T17:25:40', format='isot').mjd,
            Time('2019-06-03T17:25:56', format='isot').mjd,
            engdb_url='http://nonexistent.fake.example'
        )


def test_strict_pointing(data_file, eng_db_jw703):
    """Test failure on strict pointing"""
    with pytest.raises(ValueError):
        stp.add_wcs(data_file, siaf_path=siaf_path, tolerance=0)


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
     gs_commanded,
     fgsid,
     gs_position) = stp.get_pointing(
         Time('2019-06-03T17:25:40', format='isot').mjd,
         Time('2019-06-03T17:25:56', format='isot').mjd,
         engdb_url='http://localhost'
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
     gs_commanded,
     fgsid,
     gs_position) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, engdb_url='http://localhost')
    assert np.isclose(q, Q_EXPECTED).all()
    assert np.isclose(j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
    assert np.isclose(fsmcorr, FSMCORR_EXPECTED).all()
    assert STARTTIME <= obstime <= ENDTIME


def test_logging(eng_db_ngas, caplog):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime,
     gs_commanded,
     fgsid,
     gs_position) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, engdb_url='http://localhost')
    assert 'Determining pointing between observations times' in caplog.text
    assert 'Telemetry search tolerance' in caplog.text
    assert 'Reduction function' in caplog.text
    assert 'Querying engineering DB' in caplog.text


def test_get_pointing_list(eng_db_ngas):
    results = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, reduce_func=stp.all_pointings, engdb_url='http://localhost')
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
     gs_commanded,
     fgsid,
     gs_position) = stp.get_pointing(ZEROTIME_START.mjd, ENDTIME.mjd,
                                     reduce_func=stp.first_pointing,
                                     engdb_url='http://localhost')
    assert j2fgs_matrix.any()
    (q_desired,
     j2fgs_matrix_desired,
     fsmcorr_desired,
     obstime,
     gs_commanded,
     fgsid,
     gs_position) = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, engdb_url='http://localhost')
    assert np.array_equal(q, q_desired)
    assert np.array_equal(j2fgs_matrix, j2fgs_matrix_desired)
    assert np.array_equal(fsmcorr, fsmcorr_desired)


def test_add_wcs_default(data_file, tmp_path):
    """Handle when no pointing exists and the default is used."""
    expected_name = 'add_wcs_default.fits'

    try:
        stp.add_wcs(
            data_file, siaf_path=siaf_path, tolerance=0, allow_default=True
        )
    except ValueError:
        pass  # This is what we want for the test.
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    # Tests
    with datamodels.Level1bModel(data_file) as model:

        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        with datamodels.open(DATA_PATH / expected_name) as expected:
            for meta in METAS_EQUALITY:
                assert model[meta] == expected[meta], f'{meta} has changed'

            for meta in METAS_ISCLOSE:
                assert np.isclose(model[meta], expected[meta]), f'{meta} has changed'

            assert word_precision_check(model.meta.wcsinfo.s_region, expected.meta.wcsinfo.s_region)


def test_add_wcs_default_nosiaf(data_file_nosiaf, caplog):
    """Handle when no pointing exists and the default is used and no SIAF specified."""
    with pytest.raises(ValueError):
        stp.add_wcs(
            data_file_nosiaf, siaf_path=siaf_path, tolerance=0, allow_default=True
        )


def test_add_wcs_with_db(eng_db_ngas, data_file, tmp_path):
    """Test using the database"""
    expected_name = 'add_wcs_with_db.fits'

    stp.add_wcs(data_file, siaf_path=siaf_path, engdb_url='http://localhost')

    # Tests
    with datamodels.Level1bModel(data_file) as model:

        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        with datamodels.open(DATA_PATH / expected_name) as expected:
            for meta in METAS_EQUALITY:
                assert model[meta] == expected[meta]

            for meta in METAS_ISCLOSE:
                assert np.isclose(model[meta], expected[meta])

            assert word_precision_check(model.meta.wcsinfo.s_region, expected.meta.wcsinfo.s_region)


@pytest.mark.parametrize('fgsid', [1, 2])
def test_add_wcs_with_mast(data_file_fromsim, fgsid, tmp_path):
    """Test using the database"""
    expected_name = f'add_wcs_with_mast_fgs{fgsid}.fits'

    # See if access to MAST is available.
    try:
        engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f'Live MAST Engineering Service not available: {exception}')

    # Execute the operation.
    try:
        stp.add_wcs(data_file_fromsim, siaf_path=siaf_path, engdb_url=engdb_mast.MAST_BASE_URL,
                    fgsid=fgsid)
    except ValueError as exception:
        pytest.xfail(f'No telemetry exists. Update test to use existing telemetry. Exception: {exception}')

    # Tests
    with datamodels.Level1bModel(data_file_fromsim) as model:

        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        with datamodels.open(DATA_PATH / expected_name) as expected:
            for meta in METAS_EQUALITY:
                assert model[meta] == expected[meta], f"{meta} is not equal"

            for meta in METAS_ISCLOSE:
                assert np.isclose(model[meta], expected[meta])

            assert word_precision_check(model.meta.wcsinfo.s_region, expected.meta.wcsinfo.s_region)


def test_add_wcs_method_full_nosiafdb(eng_db_ngas, data_file, tmp_path):
    """Test using the database"""
    # Only run if `pysiaf` is installed.
    pytest.importorskip('pysiaf')

    expected_name = 'add_wcs_method_full_nosiafdb.fits'

    # Calculate
    stp.add_wcs(data_file, method=stp.Methods.OPS_TR_202111, engdb_url='http://localhost')

    # Tests
    with datamodels.Level1bModel(data_file) as model:

        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        with datamodels.open(DATA_PATH / expected_name) as expected:
            for meta in METAS_EQUALITY:
                assert model[meta] == expected[meta]

            for meta in METAS_ISCLOSE:
                assert np.isclose(model[meta], expected[meta])

            assert word_precision_check(model.meta.wcsinfo.s_region, expected.meta.wcsinfo.s_region)


def test_add_wcs_method_full_siafdb(eng_db_ngas, data_file, tmp_path):
    """Test using the database and a specified siaf db"""
    expected_name = 'add_wcs_method_full_siafdb.fits'

    # Calculate
    stp.add_wcs(data_file, siaf_path=siaf_path, method=stp.Methods.OPS_TR_202111, engdb_url='http://localhost')

    # Test
    with datamodels.Level1bModel(data_file) as model:

        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        with datamodels.open(DATA_PATH / expected_name) as expected:
            for meta in METAS_EQUALITY:
                if isinstance(model[meta], str):
                    assert model[meta] == expected[meta]
                else:
                    assert np.isclose(model[meta], expected[meta], atol=1e-13)

            for meta in METAS_ISCLOSE:
                assert np.isclose(model[meta], expected[meta])

            assert word_precision_check(model.meta.wcsinfo.s_region, expected.meta.wcsinfo.s_region)


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
        stp.update_wcs(model, siaf_path=siaf_path, allow_default=False, engdb_url='http://localhost')
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
        stp.update_wcs(model, siaf_path=siaf_path, engdb_url='http://localhost')
        assert model.meta.wcsinfo.siaf_xref_sci == 952
        assert model.meta.wcsinfo.siaf_yref_sci == 35


# ######################
# Utilities and fixtures
# ######################
def make_t_pars(detector='any', fgsid_telem=1, fgsid_user=None):
    """Setup initial Transforms Parameters

    This set was derived from the first valid group of engineering parameters for exposure
    jw00624028002_02101_00001_nrca1 retrieved from the SDP regression tests for Build 7.7.1.

    Parameters
    ==========
    fgsid_telem : [1, 2]
        The FGS reference guider to report from telemetry.

    fgsid_user : [None, 1, 2]
        The user-specified FGS to use as the reference guider.
    """
    t_pars = stp.TransformParameters()

    t_pars.detector = detector
    t_pars.fgsid = fgsid_user

    t_pars.guide_star_wcs = stp.WCSRef(ra=241.24294932221, dec=70.66165389073196, pa=None)
    t_pars.pointing = stp.Pointing(
        q=np.array([-0.20954692, -0.6177655, -0.44653177, 0.61242575]),
        j2fgs_matrix=np.array([-9.77300013e-04, 3.38988895e-03, 9.99993777e-01,
                               9.99999522e-01, 8.37175385e-09, 9.77305600e-04,
                               3.30458575e-06, 9.99994254e-01, -3.38988734e-03]),
        fsmcorr=np.array([0.00584114, -0.00432878]),
        obstime=Time(1611628160.325, format='unix'),
        gs_commanded=np.array([-22.40031242, -8.17869377]),
        gs_position=np.array([-22.4002638, -8.1786461]),
        fgsid=fgsid_telem,
    )
    t_pars.siaf = siafdb.SIAF(v2_ref=120.525464, v3_ref=-527.543132, v3yangle=-0.5627898, vparity=-1,
                              crpix1=1024.5, crpix2=1024.5, cdelt1=0.03113928, cdelt2=0.03132232,
                              vertices_idl=(-32.1682, 32.0906, 31.6586, -31.7234, -32.1683, -32.1904, 32.0823, 31.9456))
    t_pars.siaf_db = siafdb.SiafDb(siaf_path)

    return t_pars


def _calc_coarse_202111_fgsid_idfunc(value):
    """Created test IDS for calc_coarse_202111_fgsid"""
    detector, fgsid_user, fgs_expected = value
    return f'{detector}-{fgsid_user}'


def _test_methods(calc_transforms, matrix, truth_ext=''):
    """Private function to ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare

    truth_ext : str
        Arbitrary extension to add to the truth file name.
    """
    transforms, t_pars = calc_transforms

    expected_tforms = stp.Transforms.from_asdf(DATA_PATH / f'tforms_{t_pars.method}{truth_ext}.asdf')
    expected_value = getattr(expected_tforms, matrix)

    value = getattr(transforms, matrix)
    if expected_value is None:
        assert value is None
    else:
        assert np.allclose(value, expected_value)


@pytest.fixture(scope='module',
                params=(('any', None, 1), ('any', 1, 1), ('any', 2, 2),
                        ('guider1', None, 2), ('guider1', 1, 2), ('guider1', 2, 2),
                        ('guider2', None, 1), ('guider2', 1, 1), ('guider2', 2, 1)),
                ids=_calc_coarse_202111_fgsid_idfunc)
def calc_coarse_202111_fgsid(request, tmp_path_factory):
    """Calculate the transforms for COARSE_202111 with various FGS specifications
    """
    detector, fgsid_user, fgs_expected = request.param

    # Create transform parameters.
    # FGS from telemetry is set to None because, for COARSE mode,
    # telemetry is unreliable.
    t_pars = make_t_pars(detector=detector, fgsid_telem=None, fgsid_user=fgsid_user)
    t_pars.method = stp.Methods.COARSE_TR_202111
    transforms = stp.calc_transforms(t_pars)

    truth_ext = f'_{detector}-{fgsid_user}'

    # Save transforms for later examination
    transforms.write_to_asdf(tmp_path_factory.mktemp('transforms') / f'tforms_{t_pars.method}{truth_ext}.asdf')

    try:
        return transforms, t_pars, truth_ext, fgs_expected
    finally:
        t_pars.siaf_db.close()


@pytest.fixture(scope='module',
                params=[method for method in stp.Methods])
def calc_transforms(request, tmp_path_factory):
    """Calculate matrices for specified method method
    """
    t_pars = make_t_pars()

    # Set the method
    t_pars.method = request.param

    # Calculate the transforms
    transforms = stp.calc_transforms(t_pars)

    # Save transforms for later examination
    transforms.write_to_asdf(tmp_path_factory.mktemp('transforms') / f'tforms_{request.param}.asdf')

    try:
        return transforms, t_pars
    finally:
        t_pars.siaf_db.close()


@pytest.fixture
def data_file(tmp_path):
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

    file_path = tmp_path / 'file.fits'
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


@pytest.fixture
def data_file_fromsim(tmp_path):
    """Create data using times that were executed during a simulation using the OTB Simulator"""
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = Time('2022-02-02T22:24:58.942').mjd
    model.meta.exposure.end_time = Time('2022-02-02T22:26:24.836').mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.guidestar.gs_pcs_mode = 'COARSE'
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '2017-01-01'
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x_bary = -25.021
    model.meta.ephemeris.velocity_y_bary = -16.507
    model.meta.ephemeris.velocity_z_bary = -7.187

    file_path = tmp_path / 'file_fromsim.fits'
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture()
def eng_db_jw703():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_jw703_path):
        engdb = engdb_tools.ENGDB_Service(base_url='http://localhost')
        yield engdb


@pytest.fixture()
def eng_db_ngas():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_ngas_path):
        engdb = engdb_tools.ENGDB_Service(base_url='http://localhost')
        yield engdb
