"""
Unit test for Cube Build testing reading in MIRI cubepars ref file and using it
"""

import os

import numpy as np
import pytest
from astropy.io import fits
from stdatamodels.jwst.datamodels import IFUImageModel

from jwst.adaptive_trace_model import AdaptiveTraceModelStep
from jwst.adaptive_trace_model.tests import helpers
from jwst.cube_build import CubeBuildStep
from jwst.cube_build.cube_build import NoChannelsError
from jwst.cube_build.file_table import NoAssignWCSError
from jwst.datamodels import ModelContainer


@pytest.fixture(scope="module")
def miri_cube_pars(tmp_path_factory):
    """Set up the miri cube pars reference file"""

    filename = tmp_path_factory.mktemp("cube_pars")
    filename = filename / "miri_cube_pars.fits"
    hdu0 = fits.PrimaryHDU()
    hdu0.header["REFTYPE"] = "CUBEPAR"
    hdu0.header["INSTRUME"] = "MIRI"
    hdu0.header["MODELNAM"] = "FM"
    hdu0.header["DETECTOR"] = "N/A"
    hdu0.header["EXP_TYPE"] = "MIR_MRS"

    # make the first extension
    channel = np.array(["1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4"])
    subchannel = np.array(
        [
            "SHORT",
            "MEDIUM",
            "LONG",
            "SHORT",
            "MEDIUM",
            "LONG",
            "SHORT",
            "MEDIUM",
            "LONG",
            "SHORT",
            "MEDIUM",
            "LONG",
        ]
    )

    spsize = np.array([0.13, 0.13, 0.13, 0.17, 0.17, 0.17, 0.2, 0.2, 0.2, 0.35, 0.35, 0.35])
    wsamp = np.array(
        [0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.006, 0.006, 0.006]
    )

    wmin = np.array([4.89, 5.65, 6.52, 7.49, 8.65, 10.00, 11.53, 13.37, 15.44, 17.66, 20.54, 23.95])
    wmax = np.array([5.75, 6.64, 7.66, 8.78, 10.14, 11.7, 13.48, 15.63, 18.05, 20.92, 24.40, 28.45])

    col1 = fits.Column(name="CHANNEL", format="1A", array=channel)
    col2 = fits.Column(name="BAND", format="6A", array=subchannel)
    col3 = fits.Column(name="WAVEMIN", format="E", array=wmin, unit="micron")
    col4 = fits.Column(name="WAVEMAX", format="E", array=wmax, unit="micron")
    col5 = fits.Column(name="SPAXELSIZE", format="E", array=spsize, unit="arcsec")
    col6 = fits.Column(name="SPECTRALSTEP", format="D", array=wsamp, unit="micron")

    hdu1 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    hdu1.header["EXTNAME"] = "CUBEPAR"

    # make the second extension
    roispat = np.array([0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.20, 0.20, 0.20, 0.40, 0.40, 0.40])
    roispec = np.array(
        [0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.006, 0.006, 0.006]
    )

    power = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    softrad = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

    col1 = fits.Column(name="CHANNEL", format="1A", array=channel)
    col2 = fits.Column(name="BAND", format="6A", array=subchannel)
    col3 = fits.Column(name="ROISPATIAL", format="E", array=roispat, unit="arcsec")
    col4 = fits.Column(name="ROISPECTRAL", format="E", array=roispec, unit="micron")
    col5 = fits.Column(name="POWER", format="I", array=power)
    col6 = fits.Column(name="SOFTRAD", format="E", array=softrad, unit="arcsec")

    hdu2 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    hdu2.header["EXTNAME"] = "CUBEPAR_MSM"

    # make the third extension
    # Define the multiextension wavelength solution - only use a few number for testing
    finalwave = np.array([5, 10, 15, 20, 25])
    roispat = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    roispec = np.array([0.001, 0.002, 0.003, 0.004, 0.005])
    power = np.array([1, 2, 3, 4, 5])
    softrad = np.array([0.01, 0.02, 0.03, 0.04, 0.05])

    col1 = fits.Column(name="WAVELENGTH", format="D", array=finalwave, unit="micron")
    col2 = fits.Column(name="ROISPATIAL", format="E", array=roispat, unit="arcsec")
    col3 = fits.Column(name="ROISPECTRAL", format="E", array=roispec, unit="micron")
    col4 = fits.Column(name="POWER", format="I", array=power)
    col5 = fits.Column(name="SOFTRAD", format="E", array=softrad, unit="arcsec")

    hdu3 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
    hdu3.header["EXTNAME"] = "MULTICHANNEL_MSM"

    hdu = fits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdu.writeto(filename, overwrite=True)
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
