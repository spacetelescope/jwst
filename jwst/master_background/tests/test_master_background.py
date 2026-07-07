"""Unit tests for master background subtraction."""

import json
from pathlib import Path

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.cube_build import CubeBuildStep
from jwst.datamodels import ModelContainer
from jwst.extract_1d import Extract1dStep
from jwst.extract_2d import Extract2dStep
from jwst.master_background import MasterBackgroundStep
from jwst.master_background.create_master_bkg import create_background
from jwst.master_background.master_background_step import (
    copy_background_to_surf_bright,
    split_container,
)
from jwst.srctype import SourceTypeStep
from jwst.stpipe import query_step_status


@pytest.fixture(scope="module")
def user_background(tmp_path_factory):
    """Generate a user background spectrum."""
    filename = tmp_path_factory.mktemp("master_background_user_input")
    filename = filename / "user_background.fits"
    wavelength = np.linspace(0.5, 25, num=100)
    flux = np.linspace(2.0, 2.2, num=100)
    data = create_background(wavelength, flux)
    data.save(filename)
    # filename must be string not Path to validate with current stpipe version
    return str(filename)


@pytest.fixture
def nirspec_rate():
    ysize = 2048
    xsize = 2048
    shape = (ysize, xsize)
    im = datamodels.ImageModel(shape)
    im.dq = im.get_default("dq")
    im.err = im.get_default("err")
    im.var_rnoise = im.get_default("var_rnoise")
    im.var_poisson = im.get_default("var_poisson")
    im.meta.target = {"ra": 100.1237, "dec": 39.86, "source_type_apt": "EXTENDED"}
    im.meta.wcsinfo = {
        "dec_ref": 40,
        "ra_ref": 100,
        "roll_ref": 0,
        "v2_ref": -453.5134,
        "v3_ref": -373.4826,
        "v3yangle": 0.0,
        "vparity": -1,
    }
    im.meta.instrument = {
        "detector": "NRS1",
        "filter": "CLEAR",
        "grating": "PRISM",
        "name": "NIRSPEC",
        "gwa_tilt": 37.0610,
        "gwa_xtilt": 0.0001,
        "gwa_ytilt": 0.0001,
        "fixed_slit": "S200A1",
    }
    im.meta.subarray = {
        "fastaxis": 1,
        "name": "SUBS200A1",
        "slowaxis": 2,
        "xsize": 72,
        "xstart": 1,
        "ysize": 416,
        "ystart": 529,
    }
    im.meta.observation = {"program_number": "1234", "date": "2025-09-05", "time": "8:59:37"}
    im.meta.exposure = {
        "duration": 11.805952,
        "end_time": 58119.85416,
        "exposure_time": 11.776,
        "measurement_time": 11.65824,
        "frame_time": 0.11776,
        "group_time": 0.11776,
        "groupgap": 0,
        "integration_time": 11.776,
        "nframes": 1,
        "ngroups": 100,
        "nints": 1,
        "nresets_between_ints": 0,
        "nsamples": 1,
        "readpatt": "NRSRAPID",
        "sample_time": 10.0,
        "start_time": 58119.8333,
        "type": "NRS_FIXEDSLIT",
        "zero_frame": False,
    }

    yield im
    im.close()


@pytest.fixture
def nirspec_ifu_cal(nirspec_rate):
    im = datamodels.IFUImageModel(nirspec_rate)
    im.data += 1.0
    im.meta.exposure.type = "NRS_IFU"
    im.meta.subarray.name = "FULL"
    im.meta.subarray.xstart = 1
    im.meta.subarray.xsize = 2048
    im.meta.subarray.ystart = 1
    im.meta.subarray.ysize = 2048
    im = AssignWcsStep.call(im)
    yield im
    im.close()


@pytest.fixture
def nirspec_cal_pair(nirspec_rate):
    # copy the rate model to make files with different filters
    rate2 = nirspec_rate.copy()

    im1 = AssignWcsStep.call(nirspec_rate)
    im2 = AssignWcsStep.call(rate2)

    # We want to test the background subtraction, so set the
    # data to a constant.
    im1.data += 1.0
    im2.data += 1.0

    im1 = Extract2dStep.call(im1)
    im2 = Extract2dStep.call(im2)
    im1 = SourceTypeStep.call(im1)
    im2 = SourceTypeStep.call(im2)
    im1.meta.filename = "foo1_cal.fits"
    im2.meta.filename = "foo2_cal.fits"
    yield im1, im2
    im1.close()
    im2.close()


@pytest.fixture
def nirspec_ifu_pair(nirspec_ifu_cal):
    cal1 = nirspec_ifu_cal
    cal1.meta.filename = "nrs_ifu_1_cal.fits"

    s3d = CubeBuildStep.call(cal1)
    x1d = Extract1dStep.call(s3d)
    s3d.close()

    yield cal1, x1d
    cal1.close()
    x1d.close()


def _make_asn(sci_filename, x1d_filename):
    # Make a basic association with the science image
    # and the x1d as a background member
    asn = {
        "asn_type": "test",
        "asn_id": "o001",
        "asn_pool": "test",
        "products": [
            {
                "name": "product_a",
                "members": [
                    {"expname": sci_filename, "exptype": "science"},
                    {"expname": x1d_filename, "exptype": "background"},
                ],
            },
        ],
    }
    return asn


@pytest.fixture
def nirspec_asn(nirspec_cal_pair, tmp_cwd):
    """Create an association with the mock data."""
    sci, im2 = nirspec_cal_pair

    sci_filename = sci.meta.filename
    sci.save(sci_filename)

    x1d = Extract1dStep.call(im2)

    # Add a CR to the x1d
    x1d.spec[0].spec_table["SURF_BRIGHT"][10] += 10

    x1d_filename = im2.meta.filename.replace("cal", "x1d")
    x1d.meta.filename = x1d_filename
    x1d.save(x1d_filename)

    # Make a basic association with the science image
    # and the x1d as a background member
    asn = _make_asn(sci_filename, x1d_filename)

    # Save the association
    new_data = json.dumps(asn)
    asn_file = "asn.json"
    with Path(asn_file).open("w") as file:
        file.write(new_data)
    return asn_file


@pytest.fixture
def nirspec_ifu_asn(nirspec_ifu_pair, tmp_cwd):
    """Create an association with the mock data."""
    sci, x1d = nirspec_ifu_pair

    sci_filename = sci.meta.filename
    sci.save(sci_filename)

    x1d_filename = sci_filename.replace("cal", "x1d")
    x1d.meta.filename = x1d_filename
    x1d.save(x1d_filename)

    asn = _make_asn(sci_filename, x1d_filename)

    # Save the association
    new_data = json.dumps(asn)
    asn_file = "ifu_asn.json"
    with Path(asn_file).open("w") as file:
        file.write(new_data)
    return asn_file


@pytest.fixture(scope="function")
def science_image():
    """Generate science image."""
    image = datamodels.ImageModel((10, 10))
    image.dq = image.get_default("dq")
    image.err = image.get_default("err")
    image.meta.filename = "test_rate.fits"
    image.meta.instrument.name = "MIRI"
    image.meta.instrument.detector = "MIRIMAGE"
    image.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    image.meta.observation.date = "2018-01-01"
    image.meta.observation.time = "00:00:00"
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1
    image.meta.wcsinfo.v2_ref = 0
    image.meta.wcsinfo.v3_ref = 0
    image.meta.wcsinfo.v3yangle = 0
    image.meta.wcsinfo.vparity = -1
    image.meta.wcsinfo.roll_ref = 0
    image.meta.wcsinfo.ra_ref = 0
    image.meta.wcsinfo.dec_ref = 0
    image.meta.dither.x_offset = 0.0
    image.meta.dither.y_offset = 0.0
    image = AssignWcsStep.call(image)
    return image


@pytest.mark.parametrize("with_container", [True, False])
def test_master_background_userbg(tmp_cwd, user_background, science_image, with_container):
    """Verify data can run through the step with a user-supplied background."""
    if with_container:
        input_model = ModelContainer([science_image])
    else:
        input_model = science_image

    # Run with a user-supplied background and verify this is recorded in header
    output_model = MasterBackgroundStep.call(
        input_model,
        user_background=user_background,
        save_background=True,
        save_results=True,
        suffix="mbg",
    )

    if with_container:
        assert isinstance(output_model, ModelContainer)
        assert output_model is not input_model
        result = output_model[0]
    else:
        result = output_model

    assert type(science_image) is type(result)
    assert result.meta.cal_step.master_background == "COMPLETE"
    assert result.meta.background.master_background_file == "user_background.fits"

    # Input is not modified
    assert result is not science_image
    assert science_image.meta.cal_step.master_background is None

    # Expected output files are present
    assert Path("test_masterbg2d.fits").exists()
    assert Path("test_mbg.fits").exists()

    # Close models
    input_model.close()
    output_model.close()


def test_master_background_medfilt(tmp_cwd, nirspec_asn):
    """Verify data can run through the step with a median filter."""
    # Load in the example association containing a cal science member
    # and a background x1d member that has a CR in it
    container = datamodels.open(nirspec_asn)

    # Run without the median filter
    result = MasterBackgroundStep.call(
        container,
    )

    # Run with a median filter using kernel size 4. This should be rounded
    # down to 3 otherwise medfilt will throw an error.
    result_medfilt = MasterBackgroundStep.call(container, median_kernel=4)

    # Check that the background correction was run in both cases
    assert query_step_status(result, "master_background") == "COMPLETE"
    assert query_step_status(result_medfilt, "master_background") == "COMPLETE"

    # Check that there is a difference between the median filter on and off
    data_nofilt = result[0].slits[0].data
    data_medfilt = result_medfilt[0].slits[0].data
    assert not np.all(data_nofilt == data_medfilt)

    # Check that once filtered any pixels where the background correction
    # was applied have zero flux
    bkg_sub_pixels = data_medfilt[data_medfilt != 1]
    assert np.allclose(bkg_sub_pixels, 0)


@pytest.mark.parametrize("input_type", ["single", "matched_container", "mismatched_container"])
def test_master_background_logic(caplog, tmp_cwd, user_background, science_image, input_type):
    """Verify if calspec2 background step was run the master background step is skipped."""
    # the background step in calspec2 was done
    science_image.meta.cal_step.bkg_subtract = "COMPLETE"
    messages = ["Not subtracting master background"]
    if input_type == "matched_container":
        input_model = ModelContainer([science_image])
        messages.append("run again and set force_subtract")
    elif input_type == "mismatched_container":
        science_copy = science_image.copy()
        science_copy.meta.cal_step.bkg_subtract = "SKIPPED"
        input_model = ModelContainer([science_image, science_copy])
        messages.append("contains a mixture of data")
        messages.append("run again and set force_subtract")
    else:
        input_model = science_image

    # Run with a user-supplied background
    output_model = MasterBackgroundStep.call(
        input_model,
        user_background=user_background,
    )

    # Check for expected messages for the input type
    for message in messages:
        assert message in caplog.text

    if "container" in input_type:
        result = output_model[0]
    else:
        result = output_model

    assert result.meta.cal_step.master_background == "SKIPPED"
    assert type(science_image) is type(result)

    # Input is not modified
    assert result is not science_image
    assert science_image.meta.cal_step.master_background is None

    # Now force it
    output_model = MasterBackgroundStep.call(
        input_model, user_background=user_background, force_subtract=True
    )

    if "container" in input_type:
        result = output_model[0]
    else:
        result = output_model

    assert result.meta.cal_step.master_background == "COMPLETE"
    assert type(science_image) is type(result)

    # Input is still not modified
    assert result is not science_image
    assert science_image.meta.cal_step.master_background is None


def test_master_background_nirspec_ifu(tmp_cwd, nirspec_ifu_asn):
    output_model = MasterBackgroundStep.call(
        nirspec_ifu_asn, save_background=True, save_results=True, suffix="mbg"
    )

    # Output is a container
    assert isinstance(output_model, ModelContainer)
    result = output_model[0]
    assert result.meta.cal_step.master_background == "COMPLETE"

    # Expected output files are present
    basename = "nrs_ifu_1"
    expected = [
        f"{basename}_o001_masterbg1d.fits",
        f"{basename}_o001_masterbg2d.fits",
        f"{basename}_mbg.fits",
    ]
    for filename in expected:
        assert Path(filename).exists()

    # Correction should match data in most of the science regions
    # Boundaries for valid data will differ a bit from the input regions image,
    # so check the median value
    valid_regions = result.regions > 0
    assert np.isclose(np.median(result.data[valid_regions]), 0)

    # Outside those regions, all data is uncorrected
    np.testing.assert_allclose(result.data[~valid_regions], 1)

    output_model.close()


def test_master_background_invalid_input(caplog):
    input_model = datamodels.RampModel()
    result = MasterBackgroundStep.call(input_model)
    assert "cannot be handled" in caplog.text
    assert result.meta.cal_step.master_background == "FAILED"
    assert input_model.meta.cal_step.master_background is None

    input_model.close()
    result.close()


def test_master_background_no_background_data(caplog):
    input_model = ModelContainer([datamodels.ImageModel(), datamodels.ImageModel()])

    result = MasterBackgroundStep.call(input_model)
    assert "No background data found" in caplog.text
    assert result[0].meta.cal_step.master_background == "FAILED"
    assert result[1].meta.cal_step.master_background == "FAILED"

    input_model.close()
    result.close()


def test_copy_background_to_surf_bright():
    """Test the copy_background_to_surf_bright function."""
    rng = np.random.default_rng(42)
    wavelength = np.linspace(0.5, 25, num=100)
    surf_bright = np.linspace(2.0, 2.2, num=100) + (rng.random(100) - 0.5) * 0.001
    sb_error = rng.random(100) * 0.01 + 17.0  # different from berror
    background = rng.random(100) * 0.01 + 1
    berror = rng.random(100) * 0.01
    data = create_background(wavelength, surf_bright)
    data.spec[0].spec_table["sb_error"] = sb_error
    data.spec[0].spec_table["background"] = background
    data.spec[0].spec_table["bkgd_error"] = berror

    newdata = data.copy()
    copy_background_to_surf_bright(newdata)

    assert (newdata.spec[0].spec_table["surf_bright"] == background).all()
    assert (newdata.spec[0].spec_table["sb_error"] == berror).all()
    assert (newdata.spec[0].spec_table["background"] == 0).all()


def test_split_container(tmp_path):
    path1 = tmp_path / "foo.fits"
    path2 = tmp_path / "bar.fits"
    im1 = datamodels.ImageModel()
    im1.save(path1)
    im2 = datamodels.ImageModel()
    im2.save(path2)
    asn_table = {
        "asn_pool": "singleton",
        "products": [
            {
                "members": [
                    {"expname": f"{path1}", "exptype": "science"},
                    {"expname": f"{path2}", "exptype": "background"},
                ]
            }
        ],
    }
    with Path(tmp_path / "tmp_asn.json").open("w") as f:
        json.dump(asn_table, f)

    with pytest.warns(UserWarning, match="Input association file contains path information"):
        container = datamodels.open(tmp_path / "tmp_asn.json")

    sci, bkg = split_container(container)

    assert sci[0].meta.filename == "foo.fits"
    assert bkg[0].meta.filename == "bar.fits"
    assert len(sci) == 1
    assert len(bkg) == 1


def test_skip_unsupported_type(caplog):
    model = datamodels.SlitModel()
    result = MasterBackgroundStep.call(model)
    assert result.meta.cal_step.master_background == "FAILED"
    assert "SlitModel'> cannot be handled" in caplog.text

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.master_background is None


def test_skip_unsupported_without_user_bg(caplog):
    model = datamodels.ImageModel()
    result = MasterBackgroundStep.call(model)
    assert result.meta.cal_step.master_background == "FAILED"
    assert "ImageModel'> cannot be handled without user-supplied" in caplog.text

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.master_background is None
