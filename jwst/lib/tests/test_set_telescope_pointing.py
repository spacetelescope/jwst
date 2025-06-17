"""
Test suite for set_telescope_pointing
"""

import logging
import numpy as np
from pathlib import Path
import warnings

import pytest

pytest.importorskip("pysiaf")

from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.time import Time  # noqa: E402

from stdatamodels.jwst import datamodels  # noqa: E402

from jwst.lib import engdb_mast  # noqa: E402
from jwst.lib import set_telescope_pointing as stp  # noqa: E402
from jwst.lib import siafdb  # noqa: E402
from jwst.tests.helpers import word_precision_check  # noqa: E402

# Ensure that `set_telescope_pointing` logs.
stp.logger.setLevel(logging.DEBUG)
stp.logger.addHandler(logging.StreamHandler())

# Setup mock engineering service
STARTTIME = Time("2022-06-03T17:25:40", format="isot")
ENDTIME = Time("2022-06-03T17:25:56", format="isot")
ZEROTIME_START = Time("2014-01-01")
ZEROTIME_END = Time("2014-01-02")

# Header defaults
TARG_RA = 345.0
TARG_DEC = -87.0

# Get the mock databases
DATA_PATH = Path(__file__).parent / "data"

Q_EXPECTED = np.array([0.37671179, 0.70705936, -0.57895271, 0.15155541])
J2FGS_EXPECTED = np.array(
    [
        -4.93963971e-05,
        4.25990115e-03,
        9.99990925e-01,
        9.99745568e-01,
        -2.25561320e-02,
        1.45472042e-04,
        2.25565470e-02,
        9.99736502e-01,
        -4.25770310e-03,
    ]
)

J2FGS_MATRIX_EXPECTED = np.array(
    [
        -4.89526137e-05,
        4.25977338e-03,
        9.99990926e-01,
        9.99745568e-01,
        -2.25561320e-02,
        1.45025485e-04,
        2.25565451e-02,
        9.99736503e-01,
        -4.25758537e-03,
    ]
)

FSMCORR_0_EXPECTED = np.array([-0.00188433, -0.21604178])
FSMCORR_EXPECTED = np.array([0.03592864, -0.18001011])
OBSTIME_EXPECTED = Time(1654277147.967113, format="unix")

# Meta attributes for test comparisons
METAS_EQUALITY = [
    "meta.visit.engdb_pointing_quality",
    "meta.wcsinfo.wcsaxes",
    "meta.wcsinfo.ctype1",
    "meta.wcsinfo.ctype2",
    "meta.wcsinfo.cunit1",
    "meta.wcsinfo.cunit2",
    "meta.wcsinfo.vparity",
]
METAS_ISCLOSE = [
    "meta.wcsinfo.crpix1",
    "meta.wcsinfo.crpix2",
    "meta.wcsinfo.crval1",
    "meta.wcsinfo.crval2",
    "meta.wcsinfo.cdelt1",
    "meta.wcsinfo.cdelt2",
    "meta.wcsinfo.pc1_1",
    "meta.wcsinfo.pc1_2",
    "meta.wcsinfo.pc2_1",
    "meta.wcsinfo.pc2_2",
    "meta.wcsinfo.roll_ref",
    "meta.wcsinfo.v2_ref",
    "meta.wcsinfo.v3_ref",
    "meta.wcsinfo.v3yangle",
    "meta.wcsinfo.ra_ref",
    "meta.wcsinfo.dec_ref",
    "meta.pointing.ra_v1",
    "meta.pointing.dec_v1",
    "meta.pointing.pa_v3",
]


@pytest.fixture(params=[("good_model", True), ("bad_model", False), ("fits_nomodel", False)])
def file_case(request, tmp_path):
    """Generate files with different model states"""
    case, allow = request.param

    if case == "good_model":
        # Make a model that will always succeed
        model = datamodels.Level1bModel((10, 10, 10, 10))
        path = tmp_path / "level1bmodel.fits"
        model.save(path)
    elif case == "bad_model":
        # Make a model that will fail if not allowed
        model = datamodels.IFUCubeModel((10, 10, 10))
        path = tmp_path / "image.fits"
        model.save(path)
    elif case == "fits_nomodel":
        # Create just a plain anything FITS
        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        path = tmp_path / "empty.fits"
        hdul.writeto(path)
    else:
        assert False, f"Cannot produce a file for {case}"

    return path, allow


@pytest.mark.parametrize("allow_any_file", [True, False])
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
                stp.add_wcs(path, allow_any_file=allow_any_file, dry_run=True)
        else:
            # Expected error when trying to actually add the wcs.
            # The provided files do not have sufficient info to do the calculations.
            with pytest.raises(AttributeError):
                stp.add_wcs(path, allow_any_file=allow_any_file, dry_run=True)


@pytest.mark.parametrize("method", [method for method in stp.Methods])
def test_method_string(method):
    """Ensure that the value of the method is the string representation"""
    assert f"{method}" == method.value


def test_override_calc_wcs():
    """Test matrix override in the full calculation"""
    t_pars = make_t_pars()
    t_pars.method = stp.Methods.OPS_TR_202111
    wcsinfo, vinfo, _ = stp.calc_wcs(t_pars)

    override = stp.Transforms(
        m_eci2j=np.array(
            [
                [0.80583682, 0.51339893, 0.29503999],
                [-0.56953229, 0.8083677, 0.14891175],
                [-0.16204967, -0.28803339, 0.94380971],
            ]
        )
    )
    t_pars.override_transforms = override
    wcsinfo_new, vinfo_new, transforms_new = stp.calc_wcs(t_pars)

    assert vinfo_new != vinfo
    assert all(
        np.isclose(
            vinfo_new,
            stp.WCSRef(ra=32.504066294908846, dec=17.16116653015435, pa=352.2757851954435),
        )
    )


@pytest.mark.parametrize(
    "attribute, expected", [("m_eci2j", "overridden"), ("m_j2fgs1", "untouched")]
)
def test_override(attribute, expected):
    """Test overriding of Transforms attributes"""
    overrides = stp.Transforms(m_eci2j="overridden")
    to_override = stp.Transforms(m_eci2j="original", m_j2fgs1="untouched", override=overrides)

    assert getattr(to_override, attribute) == expected


def test_transform_serialize(calc_transforms, tmp_path):
    """Test serialization of Transforms"""
    transforms, t_pars = calc_transforms

    path = tmp_path / "transforms.asdf"
    transforms.write_to_asdf(path)
    from_asdf = stp.Transforms.from_asdf(path)

    assert isinstance(from_asdf, stp.Transforms)
    assert str(transforms) == str(from_asdf)


@pytest.mark.parametrize("matrix", [matrix for matrix in stp.Transforms()._fields])
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


@pytest.mark.parametrize("matrix", [matrix for matrix in stp.Transforms()._fields])
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
            Time("2015-06-15").mjd, Time("2015-06-17").mjd, engdb_url=engdb_mast.MAST_BASE_URL
        )


def test_change_engdb_url_fail():
    """Test changing the engineering database by call"""
    with pytest.raises(Exception):
        stp.get_pointing(
            Time("2019-06-03T17:25:40", format="isot").mjd,
            Time("2019-06-03T17:25:56", format="isot").mjd,
            engdb_url="http://nonexistent.fake.example",
        )


def test_strict_pointing(data_file_strict):
    """Test failure on strict pointing"""
    with pytest.raises(ValueError):
        stp.add_wcs(data_file_strict, tolerance=0)


def test_get_pointing():
    """Ensure that the averaging works."""

    (q, j2fgs_matrix, fsmcorr, obstime, gs_commanded, fgsid, gs_position) = stp.get_pointing(
        STARTTIME.mjd, ENDTIME.mjd
    )

    assert np.allclose(q, Q_EXPECTED)
    assert np.allclose(j2fgs_matrix, J2FGS_EXPECTED)
    assert np.allclose(fsmcorr, FSMCORR_EXPECTED)
    assert np.isclose(obstime.value, OBSTIME_EXPECTED.value)


def test_get_pointing_fail():
    with pytest.raises(Exception):
        (q, j2fgs_matrix, fmscorr, obstime, gs_commanded) = stp.get_pointing(47892.0, 48256.0)


def test_logging(caplog):
    (q, j2fgs_matrix, fsmcorr, obstime, gs_commanded, fgsid, gs_position) = stp.get_pointing(
        STARTTIME.mjd, ENDTIME.mjd
    )
    assert "Determining pointing between observations times" in caplog.text
    assert "Telemetry search tolerance" in caplog.text
    assert "Reduction function" in caplog.text
    assert "Querying engineering DB" in caplog.text


def test_get_pointing_list():
    results = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, reduce_func=stp.all_pointings)
    assert isinstance(results, list)
    assert len(results) > 0
    assert np.isclose(results[0].q, Q_EXPECTED).all()
    assert np.isclose(results[0].j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
    assert np.isclose(results[0].fsmcorr, FSMCORR_0_EXPECTED).all()
    assert STARTTIME <= results[0].obstime <= ENDTIME


def test_add_wcs_default(data_file, tmp_path):
    """Handle when no pointing exists and the default is used."""
    expected_name = "add_wcs_default.fits"

    try:
        stp.add_wcs(data_file, tolerance=0, allow_default=True)
    except ValueError:
        pass  # This is what we want for the test.
    except Exception as e:
        pytest.skip("Live ENGDB service is not accessible.\nException={}".format(e))

    # Tests
    with datamodels.Level1bModel(data_file) as model:
        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        with datamodels.open(DATA_PATH / expected_name) as expected:
            for meta in METAS_EQUALITY:
                assert model[meta] == expected[meta], f"{meta} has changed"

            for meta in METAS_ISCLOSE:
                assert np.isclose(model[meta], expected[meta]), f"{meta} has changed"

            assert word_precision_check(model.meta.wcsinfo.s_region, expected.meta.wcsinfo.s_region)


def test_add_wcs_default_fgsacq(tmp_path):
    """Handle when no pointing exists and the default is used."""
    with datamodels.Level1bModel(DATA_PATH / "add_wcs_default_acq1.fits") as model:
        expected_name = "add_wcs_default_acq1.fits"
        try:
            stp.update_wcs(model, tolerance=0, allow_default=True)
        except ValueError:
            pass  # This is what we want for the test.
        except Exception as e:
            pytest.skip("Live ENGDB service is not accessible.\nException={}".format(e))

        # Save for post-test comparison and update
        model.save(tmp_path / expected_name)

        # Check some keywords
        assert model.meta.wcsinfo.crpix1 == 0
        assert model.meta.wcsinfo.crpix2 == 0


def test_add_wcs_default_nosiaf(data_file_nosiaf, caplog):
    """Handle when no pointing exists and the default is used and no SIAF specified."""
    with pytest.raises(ValueError):
        stp.add_wcs(data_file_nosiaf, tolerance=0, allow_default=True)


def test_add_wcs_with_db(data_file, tmp_path):
    """Test using the database"""
    expected_name = "add_wcs_with_db.fits"

    stp.add_wcs(data_file)

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


@pytest.mark.parametrize("fgsid", [1, 2])
def test_add_wcs_with_mast(data_file_fromsim, fgsid, tmp_path):
    """Test using the database"""
    expected_name = f"add_wcs_with_mast_fgs{fgsid}.fits"

    # See if access to MAST is available.
    try:
        engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f"Live MAST Engineering Service not available: {exception}")

    # Execute the operation.
    try:
        stp.add_wcs(data_file_fromsim, engdb_url=engdb_mast.MAST_BASE_URL, fgsid=fgsid)
    except ValueError as exception:
        pytest.xfail(
            f"No telemetry exists. Update test to use existing telemetry. Exception: {exception}"
        )

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


def test_add_wcs_method_full_nosiafdb(data_file, tmp_path):
    """Test using the database"""
    expected_name = "add_wcs_method_full_nosiafdb.fits"

    # Calculate
    stp.add_wcs(data_file, method=stp.Methods.OPS_TR_202111)

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


def test_default_siaf_values(data_file_nosiaf):
    """
    Test that FITS WCS default values were set.
    """
    with datamodels.Level1bModel(data_file_nosiaf) as model:
        model.meta.exposure.start_time = STARTTIME.mjd
        model.meta.exposure.end_time = ENDTIME.mjd
        model.meta.target.ra = TARG_RA
        model.meta.target.dec = TARG_DEC
        model.meta.aperture.name = "MIRIM_TAFULL"
        model.meta.observation.date = "2017-01-01"
        model.meta.exposure.type = "MIR_IMAGE"
        stp.update_wcs(model, allow_default=False)
        assert model.meta.wcsinfo.crpix1 == 24.5
        assert model.meta.wcsinfo.crpix2 == 24.5
        assert model.meta.wcsinfo.cdelt1 == 3.0693922222222226e-05
        assert model.meta.wcsinfo.cdelt2 == 3.092664722222222e-05


def test_tsgrism_siaf_values(data_file_nosiaf):
    """
    Test that FITS WCS default values were set.
    """
    with datamodels.Level1bModel(data_file_nosiaf) as model:
        model.meta.exposure.start_time = STARTTIME.mjd
        model.meta.exposure.end_time = ENDTIME.mjd
        model.meta.aperture.name = "NRCA5_GRISM256_F444W"
        model.meta.observation.date = "2017-01-01"
        model.meta.exposure.type = "NRC_TSGRISM"
        model.meta.visit.tsovisit = True
        stp.update_wcs(model)
        assert model.meta.wcsinfo.siaf_xref_sci == 858.0
        assert model.meta.wcsinfo.siaf_yref_sci == 35


def test_mirim_tamrs_siaf_values(data_file_nosiaf):
    """
    Test that FITS WCS default values were set.
    """
    with datamodels.Level1bModel(data_file_nosiaf) as model:
        model.meta.exposure.start_time = STARTTIME.mjd
        model.meta.exposure.end_time = ENDTIME.mjd
        model.meta.aperture.name = "MIRIM_TAMRS"
        model.meta.observation.date = "2017-01-01"
        model.meta.exposure.type = "MIR_TACQ"
        stp.update_wcs(model)
        assert model.meta.wcsinfo.crpix1 == 997.5
        assert model.meta.wcsinfo.crpix2 == 993.5


def test_moving_target_update(caplog, data_file_moving_target):
    """Test moving target table updates."""
    with datamodels.open(data_file_moving_target) as model:
        stp.update_mt_kwds(model)

        expected_ra = 6.200036603575057e-05
        expected_dec = 1.7711407285091175e-10
        assert np.isclose(model.meta.wcsinfo.mt_ra, expected_ra)
        assert np.isclose(model.meta.wcsinfo.mt_dec, expected_dec)
        assert np.isclose(model.meta.target.ra, expected_ra)
        assert np.isclose(model.meta.target.dec, expected_dec)

    assert "Moving target RA and Dec updated" in caplog.text


def test_moving_target_no_mtt(caplog, data_file):
    """Test moving target table updates, no table present."""
    with datamodels.open(data_file) as model:
        stp.update_mt_kwds(model)

        # No update without table
        assert model.meta.wcsinfo.mt_ra is None
        assert model.meta.wcsinfo.mt_dec is None

    assert "Moving target position table not found" in caplog.text


def test_moving_target_tnotinrange(caplog, data_file_moving_target):
    """Test moving target table updates, time not in range."""
    with datamodels.open(data_file_moving_target) as model:
        model.meta.exposure.mid_time -= 0.2
        stp.update_mt_kwds(model)

        # No update without times in range
        assert model.meta.wcsinfo.mt_ra is None
        assert model.meta.wcsinfo.mt_dec is None

    assert "is not in the moving_target table range" in caplog.text


# ######################
# Utilities and fixtures
# ######################
def make_t_pars(detector="any", fgsid_telem=1, fgsid_user=None):
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
        j2fgs_matrix=np.array(
            [
                -9.77300013e-04,
                3.38988895e-03,
                9.99993777e-01,
                9.99999522e-01,
                8.37175385e-09,
                9.77305600e-04,
                3.30458575e-06,
                9.99994254e-01,
                -3.38988734e-03,
            ]
        ),
        fsmcorr=np.array([0.00584114, -0.00432878]),
        obstime=Time(1611628160.325, format="unix"),
        gs_commanded=np.array([-22.40031242, -8.17869377]),
        gs_position=np.array([-22.4002638, -8.1786461]),
        fgsid=fgsid_telem,
    )
    t_pars.siaf = siafdb.SIAF(
        v2_ref=120.525464,
        v3_ref=-527.543132,
        v3yangle=-0.5627898,
        vparity=-1,
        crpix1=1024.5,
        crpix2=1024.5,
        cdelt1=0.03113928,
        cdelt2=0.03132232,
        vertices_idl=(-32.1682, 32.0906, 31.6586, -31.7234, -32.1683, -32.1904, 32.0823, 31.9456),
    )
    t_pars.siaf_db = siafdb.SiafDb()

    return t_pars


def _calc_coarse_202111_fgsid_idfunc(value):
    """Created test IDS for calc_coarse_202111_fgsid"""
    detector, fgsid_user, fgs_expected = value
    return f"{detector}-{fgsid_user}"


def _test_methods(calc_transforms, matrix, truth_ext=""):
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

    expected_tforms = stp.Transforms.from_asdf(
        DATA_PATH / f"tforms_{t_pars.method}{truth_ext}.asdf"
    )
    expected_value = getattr(expected_tforms, matrix)

    value = getattr(transforms, matrix)
    if expected_value is None:
        assert value is None
    else:
        assert np.allclose(value, expected_value)


@pytest.fixture(
    scope="module",
    params=(
        ("any", None, 1),
        ("any", 1, 1),
        ("any", 2, 2),
        ("guider1", None, 2),
        ("guider1", 1, 2),
        ("guider1", 2, 2),
        ("guider2", None, 1),
        ("guider2", 1, 1),
        ("guider2", 2, 1),
    ),
    ids=_calc_coarse_202111_fgsid_idfunc,
)
def calc_coarse_202111_fgsid(request, tmp_path_factory):
    """Calculate the transforms for COARSE_202111 with various FGS specifications"""
    detector, fgsid_user, fgs_expected = request.param

    # Create transform parameters.
    # FGS from telemetry is set to None because, for COARSE mode,
    # telemetry is unreliable.
    t_pars = make_t_pars(detector=detector, fgsid_telem=None, fgsid_user=fgsid_user)
    t_pars.method = stp.Methods.COARSE_TR_202111
    transforms = stp.calc_transforms(t_pars)

    truth_ext = f"_{detector}-{fgsid_user}"

    # Save transforms for later examination
    transforms.write_to_asdf(
        tmp_path_factory.mktemp("transforms") / f"tforms_{t_pars.method}{truth_ext}.asdf"
    )

    return transforms, t_pars, truth_ext, fgs_expected


@pytest.fixture(scope="module", params=[method for method in stp.Methods])
def calc_transforms(request, tmp_path_factory):
    """Calculate matrices for specified method method"""
    t_pars = make_t_pars()

    # Set the method
    t_pars.method = request.param

    # Calculate the transforms
    transforms = stp.calc_transforms(t_pars)

    # Save transforms for later examination
    transforms.write_to_asdf(tmp_path_factory.mktemp("transforms") / f"tforms_{request.param}.asdf")

    return transforms, t_pars


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
    model.meta.observation.date = "2017-01-01"
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x = -25.021
    model.meta.ephemeris.velocity_y = -16.507
    model.meta.ephemeris.velocity_z = -7.187

    file_path = tmp_path / "file.fits"
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_strict(tmp_path):
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = STARTTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = "2017-01-01"
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x = -25.021
    model.meta.ephemeris.velocity_y = -16.507
    model.meta.ephemeris.velocity_z = -7.187

    file_path = tmp_path / "file.fits"
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_nosiaf(tmp_path):
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.aperture.name = "UNKNOWN"
    model.meta.observation.date = "2017-01-01"

    file_path = tmp_path / "fits_nosiaf.fits"
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_fromsim(tmp_path):
    """Create data using times that were executed during a simulation using the OTB Simulator"""
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = Time("2022-02-02T22:24:58.942").mjd
    model.meta.exposure.end_time = Time("2022-02-02T22:26:24.836").mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.guidestar.gs_pcs_mode = "COARSE"
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = "2017-01-01"
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x_bary = -25.021
    model.meta.ephemeris.velocity_y_bary = -16.507
    model.meta.ephemeris.velocity_z_bary = -7.187

    file_path = tmp_path / "file_fromsim.fits"
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_acq1(tmp_path):
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.aperture.name = "FGS2_FULL"
    model.meta.observation.date = "2017-01-01"
    model.meta.exposure.type = "FGS_ACQ1"
    model.meta.ephemeris.velocity_x = -25.021
    model.meta.ephemeris.velocity_y = -16.507
    model.meta.ephemeris.velocity_z = -7.187

    file_path = tmp_path / "file.fits"
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_moving_target(tmp_path):
    """Example data from simulation."""
    # Values are from simulated data file jw00634_nrcblong_mttest_uncal.fits
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = 58738.82598848102
    model.meta.exposure.end_time = 58738.82747969907
    model.meta.exposure.mid_time = 58738.82673409005
    model.meta.target.ra = 0.0
    model.meta.target.dec = 0.0
    model.meta.guidestar.gs_ra = 0.0001
    model.meta.guidestar.gs_dec = 0.0001
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = "2019-09-12"
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x = 0.00651191175424979
    model.meta.ephemeris.velocity_y = 0.160769793796114
    model.meta.ephemeris.velocity_z = 0.147663026601154

    model.meta.target.type = "MOVING"
    model.meta.moving_target = None

    times = ["2019-09-12T19:49:25.405", "2019-09-12T19:50:29.825", "2019-09-12T19:51:34.246"]
    apparent_ra = [0.0, 6.2e-5, 1.24e-4]
    apparent_dec = [-6.2e-5, 0.0, 3.0e-5]
    default = [0.0, 0.0, 0.0]
    col_names = [item["name"] for item in model.schema["properties"]["moving_target"]["datatype"]]
    mt_table = Table(
        [times, apparent_ra, apparent_dec], names=("time", "mt_apparent_RA", "mt_apparent_Dec")
    )
    for column in col_names:
        if column not in {"time", "mt_apparent_RA", "mt_apparent_Dec"}:
            mt_table.add_column(default, name=column)
    model.moving_target = mt_table.as_array()

    file_path = tmp_path / "file.fits"
    model.save(file_path)
    model.close()
    yield file_path
