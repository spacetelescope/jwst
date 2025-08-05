import logging
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels

from jwst.datamodels.utils.tso_multispec import make_tso_specmodel
from jwst.extract_1d.extract import populate_time_keywords
from jwst.tests.helpers import LogWatcher
from jwst.white_light.white_light import _determine_wavelength_range, white_light
from jwst.white_light.white_light_step import WhiteLightStep

TIME_KEYS = [
    "MJD-BEG",
    "MJD-AVG",
    "MJD-END",
    "TDB-BEG",
    "TDB-MID",
    "TDB-END",
]


@pytest.fixture(scope="module")
def make_datamodel():
    """Make data for white light tests"""
    n_spec = 5

    model = datamodels.TSOMultiSpecModel()
    model.meta.exposure.group_time = 0.15904
    model.meta.exposure.ngroups = 60
    model.meta.exposure.start_time = 58627.0
    model.meta.exposure.integration_start = 1
    model.meta.exposure.integration_end = 2
    model.meta.exposure.nints = n_spec

    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.exposure.type = "NIS_SOSS"
    model.meta.instrument.filter = "CLEAR"
    model.meta.instrument.pupil = "GR700XD"
    model.meta.observation.date = "2025-07-15"
    model.meta.observation.time = "00:00:00.000"

    # Make data arrays
    n_points = 20
    flux = np.arange(1, n_points + 1, dtype=np.float32)
    wavelength = np.linspace(11, 13, n_points, dtype=np.float32)
    f_var_poisson = np.random.rand(n_points)
    f_var_rnoise = np.random.rand(n_points)
    f_var_flat = np.random.rand(n_points)
    error = np.sqrt(f_var_poisson + f_var_rnoise + f_var_flat)

    surf_bright = np.zeros(n_points)
    sb_error = np.zeros(n_points)
    sb_var_poisson = sb_error.copy()
    sb_var_rnoise = sb_error.copy()
    sb_var_flat = sb_error.copy()
    dq = np.ones(n_points)
    background = np.zeros(n_points)
    berror = np.zeros(n_points)
    b_var_poisson = sb_error.copy()
    b_var_rnoise = sb_error.copy()
    b_var_flat = sb_error.copy()
    npixels = np.zeros(n_points)

    # This data type is used for creating an output table.
    spec_dtype = datamodels.SpecModel().spec_table.dtype

    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                f_var_poisson,
                f_var_rnoise,
                f_var_flat,
                surf_bright,
                sb_error,
                sb_var_poisson,
                sb_var_rnoise,
                sb_var_flat,
                dq,
                background,
                berror,
                b_var_poisson,
                b_var_rnoise,
                b_var_flat,
                npixels,
            ),
        ),
        dtype=spec_dtype,
    )

    spec_model = datamodels.SpecModel(spec_table=otab)
    spec_model.spectral_order = 1

    spec_list_1 = [spec_model.copy() for _ in range(n_spec)]

    # change spectral order to test multiple orders in output table
    spec_model.spectral_order = 2
    spec_list_2 = [spec_model.copy() for _ in range(n_spec)]

    times = np.linspace(0, 1, n_spec + 1)
    start_times = times[:-1]
    end_times = times[1:]
    mid_times = (start_times + end_times) / 2.0
    integrations = []
    for i in range(n_spec):
        integrations.append(
            (
                i + 1,
                start_times[i],
                mid_times[i],
                end_times[i],
                start_times[i] + 3.0,
                mid_times[i] + 3.0,
                end_times[i] + 3.0,
            )
        )

    integration_table = np.array(
        integrations,
        dtype=[
            ("integration_number", "i4"),
            ("int_start_MJD_UTC", "f8"),
            ("int_mid_MJD_UTC", "f8"),
            ("int_end_MJD_UTC", "f8"),
            ("int_start_BJD_TDB", "f8"),
            ("int_mid_BJD_TDB", "f8"),
            ("int_end_BJD_TDB", "f8"),
        ],
    )

    # make the tso spec models: one per order
    tso_spec_model_1 = make_tso_specmodel(spec_list_1)
    model.spec.append(tso_spec_model_1)
    tso_spec_model_2 = make_tso_specmodel(spec_list_2)
    model.spec.append(tso_spec_model_2)

    # update the integration times in the table
    input_model = datamodels.CubeModel((n_spec, 10, 10))
    input_model.update(model, only="PRIMARY")
    input_model.int_times = integration_table
    populate_time_keywords(input_model, model)

    # set blank times for one integration in order 1 to test warning raise
    model.spec[0].spec_table["INT_NUM"][2] = 0

    for key in TIME_KEYS:
        model.spec[0].spec_table[key][2] = np.nan

    # set one time to a different value for order 2 to test table integration
    for key in TIME_KEYS:
        model.spec[1].spec_table[key][2] += 0.01

    return model


def test_white_light(make_datamodel, monkeypatch):
    """Test white light step"""
    data = make_datamodel

    watcher = LogWatcher("1 spectra in order 1 with no mid time or duplicate mid time (20")
    monkeypatch.setattr(logging.getLogger("jwst.white_light.white_light"), "warning", watcher)
    result = white_light(data)
    watcher.assert_seen()

    n_spec = len(data.spec[0].spec_table)
    times = np.linspace(0, 1, n_spec + 1)
    start_times = times[:-1]
    end_times = times[1:]
    mid_times = (start_times + end_times) / 2.0

    # remove the mid time for the skipped integration,
    # add the adjusted midtime for order 2
    mid_times[2] = 0.51
    expected_tdb = mid_times + 3.0

    np.testing.assert_allclose(result["MJD_UTC"], mid_times)
    np.testing.assert_allclose(result["BJD_TDB"], expected_tdb, equal_nan=True)
    assert result["whitelight_flux_order_1"].shape == (len(mid_times),)
    assert result["whitelight_flux_order_2"].shape == (len(mid_times),)
    assert np.sum(np.isnan(result["whitelight_flux_order_1"])) == 1

    # check fluxes are summed appropriately
    expected_flux_per_spec = np.sum(np.arange(1, 21, dtype=np.float32))
    expected_flux = np.array(
        [
            expected_flux_per_spec,
        ]
        * len(mid_times)
    )
    expected_flux_order_1 = expected_flux.copy()
    expected_flux_order_1[-3] = np.nan  # the third was order 2 only
    expected_flux_order_2 = expected_flux.copy()

    assert_allclose(result["whitelight_flux_order_1"], expected_flux_order_1, equal_nan=True)
    assert_allclose(result["whitelight_flux_order_2"], expected_flux_order_2, equal_nan=True)


def test_white_light_multi_detector(make_datamodel):
    """Test white light step on data with multiple detectors."""
    # Set the detectors in the two spec tables to different values
    data = make_datamodel.copy()
    data.spec[0].detector = "NRS1"
    data.spec[1].detector = "NRS2"

    result = white_light(data)

    # Expected values for mid times: [0.1, 0.3, 0.5, 0.7, 0.9]
    nspec = 5

    # for NRS1, remove the 0.5 mid time for the skipped integration
    # and add a NaN at the end for padding
    expected_nrs1 = np.array([0.1, 0.3, 0.7, 0.9, np.nan])
    np.testing.assert_allclose(result["MJD_UTC_NRS1"], expected_nrs1, equal_nan=True)
    np.testing.assert_allclose(result["BJD_TDB_NRS1"], expected_nrs1 + 3.0, equal_nan=True)

    # for NRS2, adjust the 0.5 mid time to 0.51; all others are present
    expected_nrs2 = np.array([0.1, 0.3, 0.51, 0.7, 0.9])
    np.testing.assert_allclose(result["MJD_UTC_NRS2"], expected_nrs2)
    np.testing.assert_allclose(result["BJD_TDB_NRS2"], expected_nrs2 + 3.0)

    # Check for the right shape and NaN placement
    assert result["whitelight_flux_order_1_NRS1"].shape == (nspec,)
    assert result["whitelight_flux_order_2_NRS1"].shape == (nspec,)
    assert np.sum(np.isnan(result["whitelight_flux_order_1_NRS1"])) == 1
    assert np.sum(np.isnan(result["whitelight_flux_order_2_NRS1"])) == nspec

    assert result["whitelight_flux_order_1_NRS2"].shape == (nspec,)
    assert result["whitelight_flux_order_2_NRS2"].shape == (nspec,)
    assert np.sum(np.isnan(result["whitelight_flux_order_1_NRS2"])) == nspec
    assert np.sum(np.isnan(result["whitelight_flux_order_2_NRS2"])) == 0

    # check fluxes are summed appropriately
    expected_flux_per_spec = np.sum(np.arange(1, 21, dtype=np.float32))
    expected_flux = np.array(
        [
            expected_flux_per_spec,
        ]
        * nspec
    )
    expected_flux_order_1 = expected_flux.copy()
    expected_flux_order_1[-1] = np.nan  # the last entry is empty
    expected_flux_order_2 = expected_flux.copy()

    assert_allclose(result["whitelight_flux_order_1_NRS1"], expected_flux_order_1, equal_nan=True)
    assert_allclose(result["whitelight_flux_order_2_NRS2"], expected_flux_order_2, equal_nan=True)


def test_white_light_duplicate_times(monkeypatch, make_datamodel):
    """Test white light step on data with duplicate time stamps."""
    # Set all times to the same value for order 1
    data = make_datamodel.copy()
    for key in TIME_KEYS:
        data.spec[0].spec_table[key][:] = data.spec[0].spec_table[key][0]

    watcher = LogWatcher("4 spectra in order 1 with no mid time or duplicate mid time")
    monkeypatch.setattr(logging.getLogger("jwst.white_light.white_light"), "warning", watcher)
    result = white_light(data)
    watcher.assert_seen()

    # Expected times match input for the order 2 spectrum
    n_spec = len(data.spec[1].spec_table)
    mid_times = data.spec[1].spec_table["MJD-AVG"]
    expected_tdb = data.spec[1].spec_table["TDB-MID"]

    np.testing.assert_allclose(result["MJD_UTC"], mid_times)
    np.testing.assert_allclose(result["BJD_TDB"], expected_tdb, equal_nan=True)
    assert result["whitelight_flux_order_1"].shape == (len(mid_times),)
    assert result["whitelight_flux_order_2"].shape == (len(mid_times),)

    # For order 1, duplicate timestamps so only the first flux is kept
    assert np.sum(np.isnan(result["whitelight_flux_order_1"])) == n_spec - 1

    # check fluxes are summed appropriately
    expected_flux_per_spec = np.sum(np.arange(1, 21, dtype=np.float32))
    expected_flux = np.array(
        [
            expected_flux_per_spec,
        ]
        * len(mid_times)
    )
    expected_flux_order_1 = expected_flux.copy()
    expected_flux_order_1[1:] = np.nan
    expected_flux_order_2 = expected_flux.copy()

    assert_allclose(result["whitelight_flux_order_1"], expected_flux_order_1, equal_nan=True)
    assert_allclose(result["whitelight_flux_order_2"], expected_flux_order_2, equal_nan=True)


@pytest.fixture
def wavelengthrange():
    """Mock the wavelengthrange info from reference files."""
    orders = [1, 1, 2, 3]
    filters = ["CLEAR", "F277W", "CLEAR", "CLEAR"]
    wl_min = [0.93, 2.41, 0.64, 0.64]
    wl_max = [2.82, 2.82, 0.83, 0.95]
    return Table(
        data=[orders, filters, wl_min, wl_max],
        names=["order", "filter", "min_wave", "max_wave"],
        dtype=[int, str, float, float],
    )


def test_determine_wavelength_range(wavelengthrange):
    """Test that the wavelength range is determined correctly."""
    # retrieve from reference file
    wl_min, wl_max = _determine_wavelength_range(1, "F277W", waverange_table=wavelengthrange)
    assert wl_min == 2.41
    assert wl_max == 2.82

    # user-specified values override reference file values
    wl_min, wl_max = _determine_wavelength_range(
        1, "F277W", waverange_table=wavelengthrange, min_wave=2.0, max_wave=3.0
    )
    assert wl_min == 2.0
    assert wl_max == 3.0

    # use default values when no reference file is provided
    wl_min, wl_max = _determine_wavelength_range(4, "CLEAR")
    assert wl_min == -1.0
    assert wl_max == 1.0e10


def test_determine_wavelength_range_no_match(wavelengthrange):
    """Test that an error is raised if no match is found."""
    with pytest.raises(
        ValueError, match="No reference wavelength range found for order 4 and filter CLEAR"
    ):
        _determine_wavelength_range(4, "CLEAR", waverange_table=wavelengthrange)


def test_determine_wavelength_range_multiple_matches(wavelengthrange):
    """Test that an error is raised if more than one match is found."""
    wavelengthrange["order"][2] = 1
    with pytest.raises(
        ValueError, match="Multiple reference wavelength ranges found for order 1 and filter CLEAR"
    ):
        _determine_wavelength_range(1, "CLEAR", waverange_table=wavelengthrange)


def test_get_reference_wavelength_range(make_datamodel):
    """Test reading of wavelength range reference file."""
    wr = WhiteLightStep()._get_reference_wavelength_range(make_datamodel)
    assert isinstance(wr, Table)
    assert len(wr) > 0
    assert "order" in wr.columns
    assert "filter" in wr.columns
    assert "min_wave" in wr.columns
    assert "max_wave" in wr.columns


def test_get_reference_wavelength_range_other_exptype(make_datamodel):
    """Test that non-SOSS exposure types return None."""
    model = make_datamodel.copy()
    model.meta.exposure.type = "NRC_TSIMAGE"
    wr = WhiteLightStep()._get_reference_wavelength_range(model)
    assert wr is None


def test_get_reference_wavelength_range_no_file(make_datamodel, monkeypatch, log_watcher):
    """Test that missing wavelength range reference files are handled."""
    model = make_datamodel.copy()
    monkeypatch.setattr(WhiteLightStep, "get_reference_file", lambda *args: "N/A")

    watcher = log_watcher(
        "jwst.white_light.white_light_step",
        message="No wavelength range reference file found",
        level="warning",
    )
    wr = WhiteLightStep()._get_reference_wavelength_range(model)
    watcher.assert_seen()
    assert wr is None


def test_call_step(make_datamodel, tmp_cwd, log_watcher):
    """Smoke test to ensure the step at least runs without error."""
    input_copy = make_datamodel.copy()

    watcher = log_watcher(
        "jwst.white_light.white_light_step",
        message="Using wavelength range reference file",
        level="info",
    )
    result = WhiteLightStep.call(make_datamodel, save_results=True)
    watcher.assert_seen()
    assert isinstance(result, Table)
    assert Path("step_WhiteLightStep_whtlt.ecsv").exists()

    # Input is not modified
    assert result is not make_datamodel
    for input_spec, input_spec_copy in zip(make_datamodel.spec, input_copy.spec):
        np.testing.assert_allclose(
            input_spec.spec_table["FLUX"], input_spec_copy.spec_table["FLUX"]
        )
