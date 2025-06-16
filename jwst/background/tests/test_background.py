"""
Unit tests for background subtraction
"""

import json

import pytest
from numpy.testing import assert_allclose

from stdatamodels.jwst import datamodels
from jwst.background import BackgroundStep


@pytest.fixture(scope="module")
def background(tmp_path_factory):
    """Generate a background image to feed to background step"""

    filename = tmp_path_factory.mktemp("background_input")
    filename = filename / "background.fits"
    with datamodels.IFUImageModel((10, 10)) as image:
        image.data[:, :] = 10
        image.meta.instrument.name = "NIRSPEC"
        image.meta.instrument.detector = "NRS1"
        image.meta.instrument.filter = "CLEAR"
        image.meta.instrument.grating = "PRISM"
        image.meta.exposure.type = "NRS_IFU"
        image.meta.observation.date = "2019-02-27"
        image.meta.observation.time = "13:37:18.548"
        image.meta.date = "2019-02-27T13:37:18.548"

        image.meta.subarray.xstart = 1
        image.meta.subarray.ystart = 1

        image.meta.subarray.xsize = image.data.shape[-1]
        image.meta.subarray.ysize = image.data.shape[-2]

        image.meta.instrument.gwa_xtilt = 0.0001
        image.meta.instrument.gwa_ytilt = 0.0001
        image.meta.instrument.gwa_tilt = 37.0610

        image.save(filename)

    return filename


@pytest.fixture(scope="module")
def mk_asn(tmp_path_factory):
    data = {
        "asn_type": "spec2",
        "asn_rule": "Asn_Lv2Image",
        "program": "01000",
        "asn_pool": "jw010000_pool.csv",
        "products": [
            {
                "name": "jw01000005001_test_mirimage",
                "members": [
                    {
                        "expname": "jw01000005001_test_mirimage_rate.fits",
                        "exptype": "science",
                        "exposerr": "null",
                    },
                    {
                        "expname": "jw01000005001_testbg_mirimage.fits",
                        "exptype": "background",
                        "exposerr": "null",
                    },
                ],
            }
        ],
    }

    tmp_path = tmp_path_factory.mktemp("asn_input")
    asn_name = str(tmp_path / "jw010000-test_spec2_00001_asn.json")
    with open(asn_name, "w") as asn:
        json.dump(data, asn)

    # Create the files for the asn
    ratefile = str(tmp_path / "jw01000005001_test_mirimage_rate.fits")
    bgfile = str(tmp_path / "jw01000005001_testbg_mirimage.fits")
    bgsubfile = str(tmp_path / "jw010000-test_spec2_00001_asn_backgroundstep.fits")
    image_value = 10.0
    background_value = 1.0
    image = miri_rate_model((100, 100), value=image_value)
    image.save(ratefile)
    background = miri_rate_model((100, 100), value=background_value)
    background.save(bgfile)
    image.data = image.data - background.data
    image.save(bgsubfile)
    image.close()
    background.close()

    return asn_name, ratefile, bgfile, bgsubfile


@pytest.fixture(scope="function")
def science_image():
    """Generate science image"""

    image = datamodels.IFUImageModel((10, 10))
    image.data[:, :] = 100
    image.meta.instrument.name = "NIRSPEC"
    image.meta.instrument.detector = "NRS1"
    image.meta.instrument.filter = "CLEAR"
    image.meta.instrument.grating = "PRISM"
    image.meta.exposure.type = "NRS_IFU"
    image.meta.observation.date = "2019-02-27"
    image.meta.observation.time = "13:37:18.548"
    image.meta.date = "2019-02-27T13:37:18.548"
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1

    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]

    image.meta.instrument.gwa_xtilt = 0.0001
    image.meta.instrument.gwa_ytilt = 0.0001
    image.meta.instrument.gwa_tilt = 37.0610

    return image


def miri_rate_model(data_shape, value=1.0):
    """
    Generate a MIRI image subarray rate or rateints image.

    Parameters
    ----------
    data_shape : tuple of int
        Shape of the rate data. 2 values for rate, 3 for rateints.
    value : float, optional
        Value to set in the data array.

    Returns
    -------
    image : DataModel
        An open datamodel containing MIRI subarray rate or rateints
        data.
    """

    if len(data_shape) == 2:
        image = datamodels.ImageModel(data_shape)
    else:
        image = datamodels.CubeModel(data_shape)

    image.data[:, :] = value
    image.meta.instrument.name = "MIRI"
    image.meta.instrument.detector = "MIRIMAGE"
    image.meta.instrument.filter = "F2100W"
    image.meta.exposure.type = "MIR_IMAGE"
    image.meta.observation.date = "2019-02-27"
    image.meta.observation.time = "13:37:18.548"
    image.meta.date = "2019-02-27T13:37:18.548"

    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1

    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]

    return image


def test_nirspec_gwa(tmp_cwd, background, science_image):
    """Verify NIRSPEC GWA logic for in the science and background"""

    # open the background to read in the GWA values
    back_image = datamodels.open(background)
    science_image.meta.instrument.gwa_xtilt = back_image.meta.instrument.gwa_xtilt
    science_image.meta.instrument.gwa_ytilt = back_image.meta.instrument.gwa_ytilt
    science_image.meta.instrument.gwa_tilt = back_image.meta.instrument.gwa_tilt

    bkg = [background]
    # Test Run with GWA values the same - confirm it runs
    # And gives the predicted result

    result = BackgroundStep.call(science_image, bkg)

    test = science_image.data - back_image.data
    assert_allclose(result.data, test)
    assert type(result) is type(science_image)
    assert result.meta.cal_step.bkg_subtract == "COMPLETE"
    back_image.close()


def test_nirspec_gwa_xtilt(tmp_cwd, background, science_image):
    """Verify NIRSPEC GWA Xtilt must be the same in the science and background image"""

    # open the background to read in the GWA values
    back_image = datamodels.open(background)
    science_image.meta.instrument.gwa_xtilt = back_image.meta.instrument.gwa_xtilt
    science_image.meta.instrument.gwa_ytilt = back_image.meta.instrument.gwa_ytilt
    science_image.meta.instrument.gwa_tilt = back_image.meta.instrument.gwa_tilt

    bkg = [background]

    # Test change xtilt
    science_image.meta.instrument.gwa_xtilt = science_image.meta.instrument.gwa_xtilt + 0.00001

    result = BackgroundStep.call(science_image, bkg)

    assert type(result) is type(science_image)
    assert result.meta.cal_step.bkg_subtract == "SKIPPED"
    back_image.close()


def test_nirspec_gwa_ytilt(tmp_cwd, background, science_image):
    """Verify NIRSPEC GWA Ytilt must be the same in the science and background image"""

    # open the background to read in the GWA values
    back_image = datamodels.open(background)
    science_image.meta.instrument.gwa_xtilt = back_image.meta.instrument.gwa_xtilt
    science_image.meta.instrument.gwa_ytilt = back_image.meta.instrument.gwa_ytilt
    science_image.meta.instrument.gwa_tilt = back_image.meta.instrument.gwa_tilt

    bkg = [background]

    # Test different ytilt
    science_image.meta.instrument.gwa_ytilt = science_image.meta.instrument.gwa_ytilt + 0.00001

    result = BackgroundStep.call(science_image, bkg)

    assert type(result) is type(science_image)
    assert result.meta.cal_step.bkg_subtract == "SKIPPED"

    back_image.close()


@pytest.mark.parametrize(
    "data_shape,background_shape",
    [
        ((10, 10), (10, 10)),
        ((10, 10), (20, 20)),
        ((2, 10, 10), (2, 10, 10)),
        ((2, 10, 10), (2, 20, 20)),
        ((2, 10, 10), (3, 10, 10)),
        ((2, 10, 10), (3, 20, 20)),
        ((3, 10, 10), (2, 10, 10)),
        ((3, 10, 10), (2, 20, 20)),
    ],
)
def test_miri_subarray_full_overlap(data_shape, background_shape):
    image_value = 10.0
    background_value = 1.0
    image = miri_rate_model(data_shape, value=image_value)
    background = miri_rate_model(background_shape, value=background_value)

    result = BackgroundStep.call(image, [background])

    assert_allclose(result.data, image_value - background_value)
    assert type(result) is type(image)
    assert result.meta.cal_step.bkg_subtract == "COMPLETE"

    image.close()
    background.close()


@pytest.mark.parametrize(
    "data_shape,background_shape",
    [
        ((20, 20), (10, 10)),
        (
            (2, 20, 20),
            (2, 10, 10),
        ),
        (
            (3, 20, 20),
            (2, 10, 10),
        ),
        (
            (2, 20, 20),
            (3, 10, 10),
        ),
    ],
)
def test_miri_subarray_partial_overlap(data_shape, background_shape):
    image_value = 10.0
    background_value = 1.0
    image = miri_rate_model(data_shape, value=image_value)
    background = miri_rate_model(background_shape, value=background_value)

    result = BackgroundStep.call(image, [background])

    assert_allclose(
        result.data[..., : background_shape[-2], : background_shape[-1]],
        image_value - background_value,
    )
    assert_allclose(result.data[..., background_shape[-2] :, :], image_value)
    assert_allclose(result.data[..., :, background_shape[-1] :], image_value)
    assert type(result) is type(image)
    assert result.meta.cal_step.bkg_subtract == "COMPLETE"

    image.close()
    background.close()


def test_asn_input(mk_asn):
    asn_file = mk_asn[0]
    result = BackgroundStep.call(asn_file)

    bg_subtracted = mk_asn[3]
    bgs = datamodels.open(bg_subtracted)

    assert_allclose(result.data, bgs.data)
    assert result.meta.cal_step.bkg_subtract == "COMPLETE"

    result.close()
    bgs.close()


def test_bg_file_list(mk_asn):
    rate_file = mk_asn[1]
    bg_file = mk_asn[2]
    result1 = BackgroundStep.call(rate_file, bg_file)
    result2 = BackgroundStep.call(rate_file, bkg_list=bg_file)
    result3 = BackgroundStep.call(rate_file, bkg_list=[bg_file, bg_file])
    result4 = BackgroundStep.call(rate_file, bkg_list=[])
    result5 = BackgroundStep.call(rate_file)

    bg_subtracted = mk_asn[3]
    bgs = datamodels.open(bg_subtracted)

    assert_allclose(result1.data, bgs.data)
    assert_allclose(result2.data, bgs.data)
    assert_allclose(result3.data, bgs.data)

    assert result1.meta.cal_step.bkg_subtract == "COMPLETE"
    assert result2.meta.cal_step.bkg_subtract == "COMPLETE"
    assert result3.meta.cal_step.bkg_subtract == "COMPLETE"
    assert result4.meta.cal_step.bkg_subtract == "SKIPPED"
    assert result5.meta.cal_step.bkg_subtract == "SKIPPED"

    result1.close()
    result2.close()
    result3.close()
