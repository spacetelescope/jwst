"""Test output filename suffixes for the cube_build step."""

import pytest

from jwst.cube_build.ifu_cube import IFUCubeData


@pytest.fixture
def cube_instance():
    """
    Initialize IFUCubeData with the bare minimum to survive __init__.
    Use standard strings and lists
    """
    return IFUCubeData(
        input_models=[],
        output_name_base="test_base",
        output_type="band",
        linear_wave=True,
        instrument="MIRI",  # Default
        list_par1=[],
        list_par2=[],
        instrument_info={},
        master_table={},
        debug_spaxel="0 0 0",  # Must be 3 integers in a string
    )


def test_define_cubename_suffix_miri_1band(cube_instance):
    """Test MIRI suffix logic for one band and channel."""
    # 1. Set the state for a MIRI Ch1 Long cube
    cube_instance.instrument = "MIRI"
    cube_instance.list_par1 = ["1"]
    cube_instance.list_par2 = ["LONG"]

    # 2. Run the function
    suffix = cube_instance.define_cubename_suffix()

    # 3. Assert correct suffix
    assert "_ch1-long" in suffix


def test_define_cubename_suffix_miri(cube_instance):
    """Test MIRI suffix logic for multiple bands and channels."""

    cube_instance.instrument = "MIRI"
    cube_instance.list_par1 = ["1", "2", "1"]
    cube_instance.list_par2 = ["SHORT", "LONG", "MEDIUM"]

    suffix = cube_instance.define_cubename_suffix()

    # MIRI logic sorts subchannels according to increasing wavelength
    assert "_ch1-2-shortmediumlong" in suffix


def test_define_cubename_suffix_nispec(cube_instance):
    """Test NIRSpec suffix logic for multiple grating and filter"""

    cube_instance.instrument = "NIRSPEC"
    cube_instance.num_bands = 2
    cube_instance.list_par1 = ["G140M", "G140H"]
    cube_instance.list_par2 = ["F070LP", "F100LP"]

    suffix = cube_instance.define_cubename_suffix()

    assert suffix == "_g140m-f070lp-g140h-f100lp"


def test_cubename_name_nispec(cube_instance):
    """Test NIRSpec logic on stripping grating from the output_name_base"""

    cube_instance.instrument = "NIRSPEC"
    cube_instance.num_bands = 2
    cube_instance.list_par1 = ["G140M", "G140M"]
    cube_instance.list_par2 = ["F070LP", "F100LP"]
    cube_instance.output_name_base = "JP111_g140m"
    suffix = cube_instance.define_cubename_suffix()

    assert cube_instance.output_name_base == "JP111"


def test_define_cubename_suffix_nirspec(cube_instance):
    """Test NIRSpec logic for 1 grating and 1 band."""

    cube_instance.instrument = "NIRSPEC"
    cube_instance.list_par1 = ["G140M"]
    cube_instance.list_par2 = ["F100LP"]

    suffix = cube_instance.define_cubename_suffix()

    # Assert grating first then filter
    assert suffix == "_g140m-f100lp"


def test_define_cubename_suffix_internal(cube_instance):
    """Test the internal_cal coordinate system suffix."""

    cube_instance.instrument = "MIRI"
    cube_instance.list_par1 = ["4"]
    cube_instance.list_par2 = ["LONG"]
    cube_instance.coord_system = "internal_cal"

    suffix = cube_instance.define_cubename_suffix()

    assert suffix == "_ch4-long_internal"

    cube_instance.instrument = "NIRSPEC"
    cube_instance.list_par1 = ["G140M"]
    cube_instance.list_par2 = ["F100LP"]
    cube_instance.coord_system = "internal_cal"

    suffix = cube_instance.define_cubename_suffix()

    assert suffix == "_g140m-f100lp_internal"
