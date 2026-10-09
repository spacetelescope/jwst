"""Test the mod code functions."""

import os

import numpy as np
from astropy.table import Table
from stdatamodels.jwst import datamodels

import jwst.regtest.file_mod_utilities as fmod


def test_trim_tso_data(tmp_cwd):
    """Test trim_tso_data function."""
    # mock uncal tso data
    n_integrations = 10
    n_groups = 10
    rows = 20
    cols = 20
    data = np.ones([n_integrations, n_groups, rows, cols])
    table = Table()
    fake_mjd = np.arange(6000.0, 6010.0).tolist()
    table["integration_number"] = list(range(n_integrations))
    table["int_start_MJD_UTC"] = fake_mjd
    table["int_mid_MJD_UTC"] = fake_mjd
    table["int_end_MJD_UTC"] = fake_mjd
    table["int_start_BJD_TDB"] = fake_mjd
    table["int_mid_BJD_TDB"] = fake_mjd
    table["int_end_BJD_TDB"] = fake_mjd
    mock_file = tmp_cwd / "mock_file_uncal.fits"
    dm = datamodels.Level1bModel(data=data, refout=data, int_times=table.as_array())
    dm.meta.exposure.integration_start = 30
    dm.meta.exposure.integration_end = 40
    dm.save(mock_file)
    dm.close()

    os.chdir(tmp_cwd)
    expected_file = tmp_cwd / "mock_file_mod_uncal.fits"
    fmod.trim_tso_data(str(mock_file), 5, 33, 3)
    with datamodels.open(expected_file) as dm:
        assert dm.meta.filename == expected_file.name
        assert dm.meta.exposure.integration_start == 33
        assert dm.meta.exposure.integration_end == 38
        assert np.shape(dm.data) == (5, 10, 20, 20)
        assert np.shape(dm.refout) == (5, 10, 20, 20)
        assert len(dm.int_times) == 5
