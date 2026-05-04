"""
Unit test for Cube Build testing testing determining the correct suffix for the filename.
"""

import pytest

from jwst.cube_build.ifu_cube import IFUCubeData


@pytest.fixture
def cube_instance():
    """
    Initialize IFUCubeData with the bare minimum to survive __init__.
    Use standard strings and lists
    """
    return IFUCubeData(
        input_models=[],  # Empty list
        output_name_base="test_base",
        output_type="band",
        linear_wave=True,
        instrument="MIRI",  # Default
        list_par1=[],  # Empty list
        list_par2=[],  # Empty list
        instrument_info={},  # Empty dict
        master_table={},  # Empty dict
        debug_spaxel="0 0 0",  # Must be 3 integers in a string
    )


def test_define_cubename_suffix_miri_1band(monkeypatch, cube_instance):
    """Test MIRI suffix logic by monkeypatching instance attributes."""
    # 1. Set the state for a MIRI Ch1 Long cube
    monkeypatch.setattr(cube_instance, "instrument", "MIRI")
    monkeypatch.setattr(cube_instance, "list_par1", ["1"])
    monkeypatch.setattr(cube_instance, "list_par2", ["LONG"])

    # 2. Run the function
    suffix = cube_instance.define_cubename_suffix()

    # 3. Assert correct suffix
    assert "_ch1-long" in suffix


def test_define_cubename_suffix_miri(monkeypatch, cube_instance):
    """Test MIRI suffix logic by monkeypatching instance attributes."""
    # 1. Set the state for a MIRI Ch1 Short & medium, Long Ch 2
    monkeypatch.setattr(cube_instance, "instrument", "MIRI")
    monkeypatch.setattr(cube_instance, "list_par1", ["1", "2", "1"])
    monkeypatch.setattr(cube_instance, "list_par2", ["SHORT", "LONG", "MEDIUM"])

    # 2. Run the function
    suffix = cube_instance.define_cubename_suffix()

    # 3. Assert (MIRI logic sorts subchannels according to increasing wavelength)
    assert "_ch1-2-shortmediumlong" in suffix


def test_define_cubename_suffix_nirspec(monkeypatch, cube_instance):
    """Test NIRSpec logic and the base name stripping feature."""
    # 1. Setup NIRSpec state for G140M and F100LP
    monkeypatch.setattr(cube_instance, "instrument", "NIRSPEC")
    monkeypatch.setattr(cube_instance, "list_par1", ["G140M"])
    monkeypatch.setattr(cube_instance, "list_par2", ["F100LP"])
    monkeypatch.setattr(cube_instance, "num_bands", 1)

    # 2. Run the function
    suffix = cube_instance.define_cubename_suffix()

    # 3. Assertions: grating first then filter
    assert suffix == "_g140m-f100lp"


def test_define_cubename_suffix_internal(monkeypatch, cube_instance):
    """Test the internal_cal coordinate system suffix."""
    monkeypatch.setattr(cube_instance, "instrument", "MIRI")
    monkeypatch.setattr(cube_instance, "list_par1", ["1"])
    monkeypatch.setattr(cube_instance, "list_par2", ["LONG"])
    monkeypatch.setattr(cube_instance, "coord_system", "internal_cal")

    suffix = cube_instance.define_cubename_suffix()

    assert suffix == "_ch1-long_internal"
