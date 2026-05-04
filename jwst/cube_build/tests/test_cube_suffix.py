"""
Unit test for Cube Build testing testing determining the correct suffix for the filename.
"""

from unittest.mock import MagicMock

import pytest

from jwst.cube_build.ifu_cube import IFUCubeData


@pytest.fixture
def mock_ifu_cube():
    """Fixture to create an IFUCubeData instance with mocked dependencies."""
    # Mock the inputs to __init__ to avoid needing real JWST data models
    with MagicMock() as mock_models, MagicMock() as mock_inst_info:
        cube = IFUCubeData(
            input_models=mock_models,
            output_name_base="test_output",
            output_type="band",
            linear_wave=True,
            instrument="MIRI",
            list_par1=["1"],
            list_par2=["SHORT"],
            instrument_info=mock_inst_info,
            master_table={},
            debug_spaxel="1 1 1",
        )
        return cube


def test_define_cubename_suffix_miri_single(mock_ifu_cube):
    """Test MIRI suffix with a single channel/subchannel."""
    mock_ifu_cube.instrument = "MIRI"
    mock_ifu_cube.list_par1 = ["1"]
    mock_ifu_cube.list_par2 = ["SHORT"]

    suffix = mock_ifu_cube.define_cubename_suffix()
    assert suffix == "_ch1-short"


def test_define_cubename_suffix_miri_multi(mock_ifu_cube):
    """Test MIRI suffix with multiple channels and sorting logic."""
    mock_ifu_cube.instrument = "MIRI"
    mock_ifu_cube.list_par1 = ["1", "2"]
    mock_ifu_cube.list_par2 = ["SHORT", "MEDIUM"]

    suffix = mock_ifu_cube.define_cubename_suffix()
    # Note: Logic sorts subchannels descending (short -> medium -> long check)
    assert "_ch1-2" in suffix
    assert "shortmedium" in suffix  # depends on alphabetical sort


def test_define_cubename_suffix_nirspec(mock_ifu_cube):
    """Test NIRSpec suffix generation."""
    mock_ifu_cube.instrument = "NIRSPEC"
    mock_ifu_cube.list_par1 = ["G140M"]
    mock_ifu_cube.list_par2 = ["F100LP"]
    mock_ifu_cube.num_bands = 1

    suffix = mock_ifu_cube.define_cubename_suffix()
    assert suffix == "_g140m-f100lp"


def test_define_cubename_suffix_internal_cal(mock_ifu_cube):
    """Test the addition of the _internal flag."""
    mock_ifu_cube.instrument = "MIRI"
    mock_ifu_cube.list_par1 = ["1"]
    mock_ifu_cube.list_par2 = ["LONG"]
    mock_ifu_cube.coord_system = "internal_cal"

    suffix = mock_ifu_cube.define_cubename_suffix()
    assert suffix.endswith("_internal")
