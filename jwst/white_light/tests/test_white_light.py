from stdatamodels.jwst import datamodels

from jwst.white_light.white_light import white_light

import numpy as np
from numpy.testing import assert_allclose
import pytest
import logging
from jwst.tests.helpers import LogWatcher


@pytest.fixture(scope="module")
def make_datamodel():
    """Make data for white light tests"""

    model = datamodels.MultiSpecModel()
    model.meta.exposure.group_time = 0.15904
    model.meta.exposure.ngroups = 60
    model.meta.exposure.start_time = 58627.0
    model.meta.exposure.integration_start = 1
    model.meta.exposure.integration_end = 2

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

    spec_dtype = datamodels.SpecModel().spec_table.dtype  # This data type is used for creating an output table.

    otab = np.array(
        list(
            zip(
                wavelength, flux, error, f_var_poisson, f_var_rnoise, f_var_flat,
                surf_bright, sb_error, sb_var_poisson, sb_var_rnoise, sb_var_flat,
                dq, background, berror, b_var_poisson, b_var_rnoise, b_var_flat,
                npixels
            ),
        ), dtype=spec_dtype
    )

    spec_model = datamodels.SpecModel(spec_table=otab)
    spec_model.spectral_order = 1
    n_spec = 5
    [model.spec.append(spec_model.copy()) for _ in range(n_spec)]

    times = np.linspace(0, 1, n_spec+1)
    start_times = times[:-1]
    end_times = times[1:]
    mid_times = (start_times + end_times) / 2.0
    integrations = []
    for i in range(n_spec):
        integrations.append((
            i+1,
            start_times[i],
            mid_times[i],
            end_times[i],
            start_times[i] + 3.0,
            mid_times[i] + 3.0,
            end_times[i] + 3.0,
        ))

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
        ]
    )

    # also populate spec_meta with the same integration times
    for i, spec in enumerate(model.spec):
        # skip one integration to test warning raise
        if i == 2:
            continue
        # change one spectral order to test spectral ordering
        if i == 3:
            spec.spectral_order = 2
        
        spec.int_num = integration_table["integration_number"][i]
        spec.start_time_mjd = integration_table["int_start_MJD_UTC"][i]
        spec.end_time_mjd = integration_table["int_end_MJD_UTC"][i]
        spec.mid_time_mjd = integration_table["int_mid_MJD_UTC"][i]

        # remove one mid_tdb to check that code still completes and sets it to NaN
        if i != 0:
            spec.start_tdb = integration_table["int_start_BJD_TDB"][i]
            spec.end_tdb = integration_table["int_end_BJD_TDB"][i]
            spec.mid_tdb = integration_table["int_mid_BJD_TDB"][i]

    return model


def test_white_light(make_datamodel, monkeypatch):
    """Test white light step"""
    data = make_datamodel

    watcher = LogWatcher("There were 1 spectra with no mid time (20")
    monkeypatch.setattr(logging.getLogger("jwst.white_light.white_light"), "warning", watcher)
    result = white_light(data)
    watcher.assert_seen()

    n_spec = len(data.spec)
    times = np.linspace(0, 1, n_spec+1)
    start_times = times[:-1]
    end_times = times[1:]
    mid_times = (start_times + end_times) / 2.0
    mid_times = mid_times[mid_times != 0.5]  # remove the mid time for the skipped integration

    expected_tdb = mid_times + 3.0
    expected_tdb[0] = np.nan  # spec did not have mid_tdb for first integration

    np.testing.assert_allclose(result["int_mid_MJD_UTC"], mid_times)
    np.testing.assert_allclose(result["int_mid_BJD_TDB"], expected_tdb, equal_nan=True)
    assert result["whitelight_flux_order_1"].shape == (len(mid_times),)
    assert result["whitelight_flux_order_2"].shape == (len(mid_times),)
    assert np.sum(np.isnan(result["whitelight_flux_order_1"])) == 1
    assert np.sum(~np.isnan(result["whitelight_flux_order_2"])) == 1

    # check fluxes are summed appropriately
    expected_flux_per_spec = np.sum(np.arange(1, 21, dtype=np.float32))
    expected_flux = np.array([expected_flux_per_spec,]*len(mid_times))
    expected_flux_order_1 = expected_flux.copy()
    expected_flux_order_1[-2] = np.nan  # the second-to-last was order 2
    expected_flux_order_2 = expected_flux.copy()
    expected_flux_order_2[:-2] = np.nan
    expected_flux_order_2[-1] = np.nan  

    assert_allclose(result["whitelight_flux_order_1"], expected_flux_order_1, equal_nan=True)
    assert_allclose(result["whitelight_flux_order_2"], expected_flux_order_2, equal_nan=True)
