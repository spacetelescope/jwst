"""
Unit test for Cube Build testing reading in MIRI cubepars ref file and using it
"""

import os

import numpy as np
import pytest
from astropy.modeling.models import Identity, Shift
from stdatamodels.jwst.datamodels import IFUImageModel

from jwst.adaptive_trace_model import AdaptiveTraceModelStep
from jwst.adaptive_trace_model.tests import helpers
from jwst.cube_build import CubeBuildStep
from jwst.cube_build.cube_build import NoChannelsError
from jwst.cube_build.file_table import NoAssignWCSError
from jwst.cube_build.tests.helpers import make_miri_cube_pars
from jwst.datamodels import ModelContainer


@pytest.fixture(scope="module")
def miri_cube_pars(tmp_path_factory):
    filename = tmp_path_factory.mktemp("cube_pars")
    filename = filename / "miri_cube_pars.fits"
    model = make_miri_cube_pars()
    model.save(filename)
    model.close()
    return filename


@pytest.fixture(scope="function")
def miri_image_no_wcs():
    image = IFUImageModel((20, 20))
    image.data = np.random.random((20, 20))
    image.meta.instrument.name = "MIRI"
    image.meta.instrument.detector = "MIRIFUSHORT"
    image.meta.exposure.type = "MIR_MRS"
    image.meta.instrument.channel = "12"
    image.meta.instrument.band = "SHORT"
    image.meta.filename = "test_miri.fits"
    return image


@pytest.fixture(scope="module")
def nirspec_data():
    model = helpers.nirspec_ifu_model()
    model.meta.filename = "test_nirspec_cal.fits"
    return model


@pytest.fixture(scope="module")
def miri_data():
    model = helpers.miri_mrs_model()
    model.meta.filename = "test_miri_cal.fits"
    return model


@pytest.mark.parametrize("as_filename", [True, False])
def test_call_cube_build_errors(tmp_cwd, miri_cube_pars, miri_image_no_wcs, tmp_path, as_filename):
    """test defaults of step are set up and user input are defined correctly"""
    miri_image = miri_image_no_wcs
    if as_filename:
        fn = tmp_path / "miri.fits"
        miri_image.save(fn)
        step_input = fn
    else:
        step_input = miri_image

    # we do not want to run the CubeBuild through to completion because
    # the image needs to be a full image and this take too much time
    # in a unit test

    # Test NoAssignWCSError is raised
    with pytest.raises(NoAssignWCSError):
        step = CubeBuildStep()
        step.override_cubepar = miri_cube_pars
        step.channel = "3"
        step.run(step_input)

    # Test some defaults to step are setup correctly and
    # if user specifies channel it is set up correctly
    step = CubeBuildStep()
    step.override_cubepar = miri_cube_pars
    step.channel = "1"

    try:
        step.run(step_input)
    except NoAssignWCSError:
        pass

    assert step.pars_input["channel"] == ["1"]
    assert step.interpolation == "drizzle"
    assert step.weighting == "drizzle"
    assert step.coord_system == "skyalign"

    # Set Assign WCS has been run but the user input to channels is wrong
    miri_image.meta.cal_step.assign_wcs = "COMPLETE"
    # save file with modifications
    if as_filename:
        miri_image.save(step_input)
    with pytest.raises(NoChannelsError):
        step = CubeBuildStep()
        step.override_cubepar = miri_cube_pars
        step.channel = "3"
        step.run(step_input)


@pytest.mark.parametrize("as_filename", [True, False])
@pytest.mark.parametrize("coord_system", ["internal_cal", "skyalign"])
def test_call_cube_build_nirspec(tmp_cwd, nirspec_data, tmp_path, as_filename, coord_system):
    # Add a NaN in the error array, unmatched in data, to
    # check that the input is not modified by the match_nans_and_flags
    # call in the beginning of the step
    nirspec_data.err[100, 100] = np.nan
    if as_filename:
        fn = tmp_path / "test_nirspec_cal.fits"
        nirspec_data.save(fn)
        step_input = fn
    else:
        step_input = nirspec_data.copy()
    step = CubeBuildStep()
    step.coord_system = coord_system
    step.save_results = True
    result = step.run(step_input)

    assert isinstance(result, ModelContainer)
    assert len(result) == 1
    model = result[0]
    assert model.meta.cal_step.cube_build == "COMPLETE"
    if coord_system == "internal_cal":
        assert model.meta.filename == "test_nirspec_prism-clear_internal_s3d.fits"
    else:
        assert model.meta.filename == "test_nirspec_prism-clear_s3d.fits"
    assert os.path.isfile(model.meta.filename)

    # make sure input is not modified
    assert result is not step_input
    assert result[0] is not step_input
    if not as_filename:
        np.testing.assert_allclose(step_input.data, nirspec_data.data)
        assert step_input.meta.cal_step.cube_build is None


def test_missing_cubepars(nirspec_data):
    with pytest.raises(ValueError, match="cubepar reference file is required"):
        CubeBuildStep.call(nirspec_data, override_cubepar="N/A")


def test_invalid_coord_sys(miri_image_no_wcs, miri_cube_pars):
    step = CubeBuildStep()
    step.override_cubepar = miri_cube_pars
    step.coord_system = "internal_cal"
    miri_image_no_wcs.meta.cal_step.assign_wcs = "COMPLETE"
    with pytest.raises(ValueError, match="coordinate system is not supported for MIRI"):
        step.run(miri_image_no_wcs)


@pytest.mark.parametrize("oversample", [1, 2])
@pytest.mark.parametrize(
    "dataset,output_shape", [("nirspec_data", (941, 41, 37)), ("miri_data", (850, 41, 41))]
)
def test_cube_build_oversampled(request, oversample, dataset, output_shape):
    model = request.getfixturevalue(dataset)

    # Oversample the input data, or just attach a profile
    shape = model.data.shape
    oversampled = AdaptiveTraceModelStep.call(model, oversample=oversample)
    if dataset.startswith("miri"):
        assert oversampled.data.shape == (shape[0], shape[1] * oversample)
    else:
        assert oversampled.data.shape == (shape[0] * oversample, shape[1])
    if oversample != 1:
        assert oversampled.hasattr("regions")

    # In either case, cube build succeeds and output size is the about the same
    cube = CubeBuildStep.call(oversampled)[0]
    assert cube.meta.cal_step.cube_build == "COMPLETE"
    np.testing.assert_allclose(cube.data.shape, output_shape, atol=2)

    # Input data is flat, so output data should have the same mean value
    np.testing.assert_allclose(np.nanmean(cube.data), np.nanmean(model.data))


@pytest.mark.parametrize("shift_ra", [True, False])
def test_output_crval1_positive(miri_data, shift_ra):
    model = miri_data.copy()

    expected_ra = 164
    if shift_ra:
        # Shift RA by 100 degrees to make it > 180
        model.meta.wcs.pipeline[-2].transform |= Shift(100) & Identity(2)
        model.meta.wcsinfo.ra_ref += 100
        model.meta.wcsinfo.s_region = model.meta.wcsinfo.s_region.replace("164", "264")
        expected_ra += 100

    result = CubeBuildStep.call(model)

    # Output RA should be between 0 and 360
    for cube in result:
        assert np.isclose(cube.meta.wcsinfo.crval1, expected_ra, atol=1)
