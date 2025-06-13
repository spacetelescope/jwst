"""
Unit test for Cube Build testing setting up configuration
"""

import pytest
import math
import asdf
from stdatamodels.jwst import datamodels
import astropy.units as u
from gwcs import WCS
import numpy as np
from jwst.cube_build import CubeBuildStep
from jwst.cube_build import cube_build
from jwst.cube_build import ifu_cube
from jwst.cube_build import file_table
from jwst.cube_build import instrument_defaults
from jwst import assign_wcs


@pytest.fixture(scope="module")
def offset_file(tmp_path_factory):
    """Generate a offset file"""

    filename = tmp_path_factory.mktemp("offset")
    filename = filename / "offset.asdf"

    testfile = ["test1.fits", "test2.fits"]
    raoffset = [0.0, 0.1]
    decoffset = [0.0, 0.15]
    tree = {
        "units": str(u.arcsec),
        "filename": testfile,
        "raoffset": raoffset,
        "decoffset": decoffset,
    }
    af = asdf.AsdfFile(tree)
    af.write_to(filename)
    af.close()
    return filename


@pytest.fixture(scope="module")
def offset_file_arcmin(tmp_path_factory):
    """Generate a offset file with units = arcmin"""

    filename = tmp_path_factory.mktemp("offset")
    filename = filename / "offset_arcmin.asdf"

    testfile = ["test1.fits", "test2.fits"]
    raoffset = [0.0, 0.1]
    decoffset = [0.0, 0.15]
    tree = {
        "units": str(u.arcmin),
        "filename": testfile,
        "raoffset": raoffset,
        "decoffset": decoffset,
    }
    af = asdf.AsdfFile(tree)
    af.write_to(filename)
    return filename


@pytest.fixture(scope="function")
def miri_ifushort_short_2files():
    """Generate input model with 2 IFU images"""

    observation = {"date": "2019-01-01", "time": "17:00:00"}

    subarray = {
        "fastaxis": 1,
        "name": "FULL",
        "slowaxis": 2,
        "xsize": 1032,
        "xstart": 1,
        "ysize": 1024,
        "ystart": 1,
    }

    wcsinfo = {
        "dec_ref": 39.05036271706514,
        "ra_ref": 339.8149235604264,
        "roll_ref": 217.25027556008598,
        "v2_ref": -503.378,
        "v3_ref": -318.9992,
        "v3yangle": 0.0,
        "vparity": -1,
        "s_region": "POLYGON ICRS  339.813915797 39.049575409 339.816080118 39.049575409 "
        + "339.816080118 39.051260090 339.813915797 39.051260090",
        "spectral_region": ([4.889451133245338, 8.781164838427532]),
    }

    mirifushort_short = {
        "detector": "MIRIFUSHORT",
        "channel": "12",
        "band": "SHORT",
        "name": "MIRI",
    }

    input_model1 = datamodels.IFUImageModel()
    input_model1.meta.exposure.type = "MIR_MRS"
    input_model1.meta.wcsinfo._instance.update(wcsinfo)
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.observation._instance.update(observation)
    input_model1.meta.subarray._instance.update(subarray)
    input_model1.meta.cal_step.assign_wcs = "COMPLETE"
    input_model1.meta.filename = "test1.fits"
    input_model1.data = np.random.random((1024, 1032))

    input_model2 = datamodels.IFUImageModel()
    input_model2.meta.exposure.type = "MIR_MRS"
    input_model2.meta.wcsinfo._instance.update(wcsinfo)
    input_model2.meta.instrument._instance.update(mirifushort_short)
    input_model2.meta.observation._instance.update(observation)
    input_model2.meta.subarray._instance.update(subarray)
    input_model2.meta.cal_step.assign_wcs = "COMPLETE"
    input_model2.meta.filename = "test2.fits"
    input_model2.data = np.random.random((1024, 1032))

    input_models = []

    step = assign_wcs.assign_wcs_step.AssignWcsStep()
    refs = {}
    for reftype in assign_wcs.assign_wcs_step.AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(input_model1, reftype)
    pipe = assign_wcs.miri.create_pipeline(input_model1, refs)
    input_model1.meta.wcs = WCS(pipe)

    for reftype in assign_wcs.assign_wcs_step.AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(input_model2, reftype)
    pipe = assign_wcs.miri.create_pipeline(input_model2, refs)
    input_model2.meta.wcs = WCS(pipe)

    input_models.append(input_model1)
    input_models.append(input_model2)
    return input_models


def test_offset_file_config(tmp_cwd, miri_ifushort_short_2files, offset_file):
    """Test validation of the offset configuration"""

    # first test that it is a valid asdf file and has what is needed
    step = CubeBuildStep()
    step.input_models = miri_ifushort_short_2files

    step.offset_file = offset_file
    offsets = step.check_offset_file()
    assert isinstance(offsets, dict)


def test2_offset_file_config(tmp_cwd, miri_ifushort_short_2files, offset_file):
    """Test validation of the offset configuration"""

    # Test changing one of the filenames so it is not in the list given
    # in the offset_file
    step = CubeBuildStep()
    step.input_models = miri_ifushort_short_2files

    miri_ifushort_short_2files[0].meta.filename = "test3.fits"
    step.offset_file = offset_file

    with pytest.raises(ValueError):
        step.check_offset_file()


def test_offset_file_units(tmp_cwd, miri_ifushort_short_2files, offset_file_arcmin):
    """Test offsets are not used when units are arc minutes"""

    # test is the if the user set the units to arcmins
    step = CubeBuildStep()
    step.input_models = miri_ifushort_short_2files

    step.offset_file = offset_file_arcmin
    with pytest.raises(Exception):
        step.check_offset_file()


def test_read_offset_file(miri_ifushort_short_2files, offset_file):
    """Test offset file has been read in correctly"""

    step = CubeBuildStep()
    step.input_models = miri_ifushort_short_2files
    step.offset_file = offset_file
    offsets = step.check_offset_file()
    # Test that the offset file is read in and is a dictionary
    assert isinstance(offsets, dict)

    pars_input = {}
    pars_input["channel"] = []
    pars_input["subchannel"] = []
    pars_input["filter"] = []
    pars_input["grating"] = []
    weighting = "drizzle"
    output_type = "multi"
    single = False
    par_filename = "None"

    # set up pars needed for CubeData class
    pars = {
        "channel": pars_input["channel"],
        "subchannel": pars_input["subchannel"],
        "grating": pars_input["grating"],
        "filter": pars_input["filter"],
        "weighting": weighting,
        "single": single,
        "output_type": output_type,
        "offset_file": offset_file,
    }

    cubeinfo = cube_build.CubeData(miri_ifushort_short_2files, par_filename, **pars)

    master_table = file_table.FileTable()
    this_instrument = master_table.set_file_table(cubeinfo.input_models)

    cubeinfo.instrument = this_instrument
    cubeinfo.determine_band_coverage(master_table)
    num_cubes, cube_pars = cubeinfo.number_cubes()
    # test with output_type = mulit we get 1 cube
    # test that cube info sets up the correct channels and band for data
    assert num_cubes == 1
    assert cube_pars["1"]["par1"] == ["1", "2"]
    assert cube_pars["1"]["par2"] == ["short", "short"]

    wave_min = 4.88
    wave_max = 8.78

    # set up par for IFUCubeData CLASS
    pars_cube = {
        "scalexy": 0.0,
        "scalew": 0.0,
        "interpolation": "drizzle",
        "weighting": "drizzle",
        "weight_power": None,
        "coord_system": "skyalign",
        "ra_center": None,
        "dec_center": None,
        "cube_pa": None,
        "nspax_x": None,
        "nspax_y": None,
        "rois": None,
        "riow": None,
        "wavemin": wave_min,
        "wavemax": wave_max,
        "skip_dqflagging": False,
        "offsets": offsets,
        "debug_spaxel": "0 0 0",
    }

    pipeline = 3
    list_par1 = ["1", "2"]
    list_par2 = ["short", "short"]
    output_name_base = "TEMP"

    instrument_info = instrument_defaults.InstrumentInfo()

    thiscube = ifu_cube.IFUCubeData(
        pipeline,
        miri_ifushort_short_2files,
        output_name_base,
        output_type,
        this_instrument,
        list_par1,
        list_par2,
        instrument_info,
        master_table,
        **pars_cube,
    )

    thiscube.linear_wavelength = True
    thiscube.spatial_size = 0.13
    thiscube.spectral_size = 0.001
    thiscube.setup_ifucube_wcs()

    # test the offset file was read in correctly
    filename = ["test1.fits", "test2.fits"]
    raoffset = [0.0, 0.1]
    decoffset = [0.0, 0.15]

    ravalues = thiscube.offsets["raoffset"]
    decvalues = thiscube.offsets["decoffset"]

    assert thiscube.offsets["filename"] == filename

    assert math.isclose(ravalues[0].value, raoffset[0], abs_tol=0.0001)
    assert math.isclose(ravalues[1].value, raoffset[1], abs_tol=0.0001)
    assert math.isclose(decvalues[0].value, decoffset[0], abs_tol=0.0001)
    assert math.isclose(decvalues[1].value, decoffset[1], abs_tol=0.0001)
