import logging

import numpy as np
from numpy.testing import assert_allclose
import pytest
from stdatamodels.jwst import datamodels

from jwst.datamodels.utils.tso_multispec import make_tso_specmodel
from jwst.extract_1d.extract import populate_time_keywords
from jwst.tests.helpers import LogWatcher
from jwst.white_light.white_light import white_light


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
    time_keys = [
        "MJD-BEG",
        "MJD-AVG",
        "MJD-END",
        "TDB-BEG",
        "TDB-MID",
        "TDB-END",
    ]
    for key in time_keys:
        model.spec[0].spec_table[key][2] = np.nan

    # set one time to a different value for order 2 to test table integration
    for key in time_keys:
        model.spec[1].spec_table[key][2] += 0.01

    return model


def test_white_light(make_datamodel, monkeypatch):
    """Test white light step"""
    data = make_datamodel

    watcher = LogWatcher("There were 1 spectra in order 1 with no mid time (20")
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
    """Test white light step"""
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
