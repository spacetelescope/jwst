"""Unit tests for master background NIRSpec corrections."""

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.master_background.nirspec_utils import (
    correct_nrs_fs_bkg,
    correct_nrs_ifu_bkg,
    is_background_msa_slit,
)

FS_CORRECTION_ARRAYS = (
    "flatfield_point",
    "flatfield_uniform",
    "pathloss_point",
    "pathloss_uniform",
    "photom_point",
    "photom_uniform",
)


@pytest.mark.parametrize("missing", ["pathloss_point", "pathloss_uniform"])
def test_ifu_pathloss_existence(missing):
    """Test the case where the input is missing a pathloss array."""
    input_data = datamodels.IFUImageModel((10, 10))
    if missing == "pathloss_point":
        input_data.pathloss_uniform = input_data.get_default("pathloss_uniform")
    else:
        input_data.pathloss_point = input_data.get_default("pathloss_point")
    input_copy = input_data.copy()
    result = correct_nrs_ifu_bkg(input_data)

    # If either array is missing, output is the same as input
    assert np.all(result.data == input_copy.data)


def test_ifu_correction():
    """Test application of IFU corrections."""
    data = np.ones((5, 5))
    pl_ps = 2.1 * data
    pl_un = data / 1.9
    input_data = datamodels.IFUImageModel(data=data, pathloss_point=pl_ps, pathloss_uniform=pl_un)

    corrected = input_data.data * pl_un / pl_ps
    result = correct_nrs_ifu_bkg(input_data)

    assert np.allclose(corrected, result.data, rtol=1.0e-7)


def test_fs_correction():
    """Test application of FS corrections."""
    data = np.ones((5, 5))
    ff_ps = 1.5 * data
    ff_un = data / 1.2
    pl_ps = 2 * data
    pl_un = data / 2.1
    ph_ps = 1.1 * data
    ph_un = 1.23 * data
    input_data = datamodels.SlitModel(
        data=data,
        flatfield_point=ff_ps,
        flatfield_uniform=ff_un,
        pathloss_point=pl_ps,
        pathloss_uniform=pl_un,
        photom_point=ph_ps,
        photom_uniform=ph_un,
    )

    corrected = input_data.data * (ff_un / ff_ps) * (pl_un / pl_ps) * (ph_ps / ph_un)
    result = correct_nrs_fs_bkg(input_data)

    assert np.allclose(corrected, result.data, rtol=1.0e-7)


@pytest.mark.parametrize("missing", FS_CORRECTION_ARRAYS)
def test_fs_correction_missing_array(caplog, missing):
    """Test FS correction with missing arrays."""
    data = np.ones((5, 5))
    input_data = datamodels.SlitModel(data=data)
    for array_name in FS_CORRECTION_ARRAYS:
        if array_name != missing:
            setattr(input_data, array_name, data.copy())

    # No correction if any array is missing
    result = correct_nrs_fs_bkg(input_data)
    assert np.allclose(result.data, 1.0)
    assert f"{missing} array not found"


@pytest.mark.parametrize(
    "name,status",
    [("BKG101", True), ("bkg101", True), ("background_101", True), ("101", False), (None, False)],
)
def test_is_background(name, status):
    """Test check for background slit."""
    slit = datamodels.SlitModel()
    if name is not None:
        slit.source_name = name
    assert is_background_msa_slit(slit) == status
