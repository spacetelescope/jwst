"""
Unit test for mrs_imatch testing setting up configuration.
"""

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.mrs_imatch.mrs_imatch_step import MRSIMatchStep, _get_2d_pixgrid, _find_channel_bkg_index

mirifushort_short = {"detector": "MIRIFUSHORT", "channel": "12", "band": "SHORT", "name": "MIRI"}

mirifulong_long = {"detector": "MIRIFULONG", "channel": "34", "band": "LONG", "name": "MIRI"}


@pytest.fixture(scope="function")
def miri_dither_ch12():
    """Generate 4 dithered channel 12 data."""

    input_model1 = datamodels.IFUImageModel((20, 20))
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.cal_step.assign_wcs = "COMPLETE"

    input_model2 = datamodels.IFUImageModel((20, 20))
    input_model2.meta.instrument._instance.update(mirifushort_short)
    input_model2.meta.cal_step.assign_wcs = "COMPLETE"

    input_model3 = datamodels.IFUImageModel((20, 20))
    input_model3.meta.instrument._instance.update(mirifushort_short)
    input_model3.meta.cal_step.assign_wcs = "COMPLETE"

    input_model4 = datamodels.IFUImageModel((20, 20))
    input_model4.meta.instrument._instance.update(mirifushort_short)
    input_model4.meta.cal_step.assign_wcs = "COMPLETE"

    # stuff in model container
    input_models = []
    input_models.append(input_model1)
    input_models.append(input_model2)
    input_models.append(input_model3)
    input_models.append(input_model4)

    return input_models


def test_imatch_background_subtracted(tmp_cwd, miri_dither_ch12):
    """Test if data is already background subtracted - raise error."""

    all_models = ModelContainer(miri_dither_ch12)
    # modify the data set background subtracted
    new_container = []
    for m in all_models:
        m.meta.background.subtracted = True
        new_container.append(m)

    # test if background subtracted - raise error
    with pytest.raises(ValueError):
        step = MRSIMatchStep()
        step.call(new_container, skip=False)


def test_imatch_default_run():
    """Test mrs_imatch test is skipped by default"""

    # test if default running results in skipping step
    step = MRSIMatchStep()
    assert step.skip


def test_imatch_background_reset(tmp_cwd, miri_dither_ch12):
    """Test if background polynomial is already determined - reset it"""

    all_models = ModelContainer(miri_dither_ch12)

    # added a background and test is reset background
    # removes the background
    new_container = []
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "2"
    for m in all_models:
        m.meta.background.polynomial_info.append(
            {
                "degree": degree,
                "refpoint": center,
                "coefficients": poly.ravel().tolist(),
                "channel": channel,
            }
        )
        new_container.append(m)

    # test if reset background - removes background
    step = MRSIMatchStep()
    step._reset_background(new_container)

    for i in range(len(new_container)):
        m = new_container[i]
        assert len(m.meta.background.polynomial_info) == 0


def test_find_channel_index(tmp_cwd, miri_dither_ch12):
    """Test if correct channel index is returned"""

    # channel 1 - model only has 1 background polynomial
    input_model12 = datamodels.IFUImageModel((20, 20))
    input_model12.meta.instrument._instance.update(mirifushort_short)
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "1"
    input_model12.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    assert _find_channel_bkg_index(input_model12, "1") == 0

    # channel 2  background only has 1 background polynomial
    input_model12 = datamodels.IFUImageModel((20, 20))
    input_model12.meta.instrument._instance.update(mirifushort_short)
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "2"
    input_model12.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    assert _find_channel_bkg_index(input_model12, "2") == 0

    # add background polynomial for channel 1 and test index
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "1"
    input_model12.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    assert _find_channel_bkg_index(input_model12, "1") == 1

    # set up a new models only channel 3
    # channel 3
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "3"
    input_model34.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    assert _find_channel_bkg_index(input_model34, "3") == 0

    # set up a new model only have channel 4
    # channel 4
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "4"
    input_model34.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    assert _find_channel_bkg_index(input_model34, "4") == 0

    # set up a new model with both channel 3 and 4 background polynomials
    # channel 3 & 4
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    degree = (
        1,
        1,
        1,
    )
    center = (
        5,
        5,
        5,
    )
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = "4"
    input_model34.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    channel = "3"
    input_model34.meta.background.polynomial_info.append(
        {
            "degree": degree,
            "refpoint": center,
            "coefficients": poly.ravel().tolist(),
            "channel": channel,
        }
    )

    assert _find_channel_bkg_index(input_model34, "4") == 0
    assert _find_channel_bkg_index(input_model34, "3") == 1


@pytest.mark.parametrize("channel", ["1", "2", "3", "4"])
def test_get_2d_pixgrid(channel):
    """Test if x,y grid for channel is formed correctly."""
    # channels 1 and 4: left side of detector
    # test shape of x (1/2 detector) and min,max (0,511) + 4 reference pixels
    #
    # channels 2 and 3; right side of detector
    # test shape of x (1/2 detector) and min,max 516, 1027 (no reference pixels)
    x, y = _get_2d_pixgrid(channel)
    assert x.shape == (1024, 512)
    if channel in ("1", "4"):
        assert np.amin(x) == 4
        assert np.amax(x) == 515
    else:
        assert np.amin(x) == 516
        assert np.amax(x) == 1027
