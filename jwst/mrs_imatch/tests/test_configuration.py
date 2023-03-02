"""
Unit test for mrs_imatch testing setting up configuration
"""

import pytest

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.mrs_imatch.mrs_imatch_step import MRSIMatchStep
from jwst.mrs_imatch.mrs_imatch_step import _get_2d_pixgrid
from jwst.mrs_imatch.mrs_imatch_step import _find_channel_bkg_index
import numpy as np

mirifushort_short = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'SHORT',
    'name': 'MIRI'
}

mirifulong_long = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'LONG',
    'name': 'MIRI'
}


@pytest.fixture(scope='function')
def miri_dither_ch12():
    """ Generate 4 dithered channel 12 data  """

    input_model1 = datamodels.IFUImageModel((20, 20))
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model2 = datamodels.IFUImageModel((20, 20))
    input_model2.meta.instrument._instance.update(mirifushort_short)
    input_model2.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model3 = datamodels.IFUImageModel((20, 20))
    input_model3.meta.instrument._instance.update(mirifushort_short)
    input_model3.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model4 = datamodels.IFUImageModel((20, 20))
    input_model4.meta.instrument._instance.update(mirifushort_short)
    input_model4.meta.cal_step.assign_wcs = 'COMPLETE'

    # stuff in model container
    input_models = []
    input_models.append(input_model1)
    input_models.append(input_model2)
    input_models.append(input_model3)
    input_models.append(input_model4)

    return input_models


def test_imatch_background_subtracted(_jail, miri_dither_ch12):
    """ Test if data is already background subtracted - raise error"""

    all_models = ModelContainer(miri_dither_ch12)
    # modify the data set background subtracted
    new_container = []
    for m in all_models:
        m.meta.background.subtracted = True
        new_container.append(m)

    # test if background subtracted - raise error
    with pytest.raises(ValueError):
        step = MRSIMatchStep()
        step.run(new_container)


def test_imatch_background_reset(_jail, miri_dither_ch12):
    """ Test if background polynomial is already determined - reset it"""

    all_models = ModelContainer(miri_dither_ch12)

    # added a background and test is reset background
    # removes the background
    new_container = []
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '2'
    for m in all_models:
        m.meta.background.polynomial_info.append(
            {
                'degree': degree,
                'refpoint': center,
                'coefficients': poly.ravel().tolist(),
                'channel': channel
            }
        )
        new_container.append(m)

    # test if reset background - removes background
    step = MRSIMatchStep()
    step._reset_background(new_container)

    for i in range(len(new_container)):
        m = new_container[i]
        test = len(m.meta.background.polynomial_info)
        assert test == 0


def test_find_channel_index(_jail, miri_dither_ch12):
    """ Test if correct channel index is returned """

    # channel 1 - model only has 1 background polynomial
    input_model12 = datamodels.IFUImageModel((20, 20))
    input_model12.meta.instrument._instance.update(mirifushort_short)
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '1'
    input_model12.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    index = _find_channel_bkg_index(input_model12, '1')
    assert index == 0

    # channel 2  background only has 1 background polynomial
    input_model12 = datamodels.IFUImageModel((20, 20))
    input_model12.meta.instrument._instance.update(mirifushort_short)
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '2'
    input_model12.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    index = _find_channel_bkg_index(input_model12, '2')
    assert index == 0

    # add background polynomial for channel 1 and test index
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '1'
    input_model12.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    index = _find_channel_bkg_index(input_model12, '1')
    assert index == 1

    # set up a new models only channel 3
    # channel 3
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '3'
    input_model34.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    index = _find_channel_bkg_index(input_model34, '3')
    assert index == 0

    # set up a new model only have channel 4
    # channel 4
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '4'
    input_model34.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    index = _find_channel_bkg_index(input_model34, '4')
    assert index == 0

    # set up a new model with both channel 3 and 4 background polynomials
    # channel 3 & 4
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    degree = (1, 1, 1,)
    center = (5, 5, 5,)
    poly = np.ndarray(9)
    poly[:] = 1.3
    channel = '4'
    input_model34.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    channel = '3'
    input_model34.meta.background.polynomial_info.append(
        {
            'degree': degree,
            'refpoint': center,
            'coefficients': poly.ravel().tolist(),
            'channel': channel
        }
    )

    index = _find_channel_bkg_index(input_model34, '4')
    assert index == 0

    index = _find_channel_bkg_index(input_model34, '3')
    assert index == 1


def test_get_2d_pixgrid():
    """ Test if x,y grid for channel is formed correctly  """

    # test if channel 1 then left side of detector
    # test shape of x (1/2 detector) and min,max (0,511) + 4 reference pixels
    input_model12 = datamodels.IFUImageModel((20, 20))
    input_model12.meta.instrument._instance.update(mirifushort_short)
    x, y = _get_2d_pixgrid(input_model12, '1')

    xsize0 = x.shape[0]
    xsize1 = x.shape[1]
    assert xsize0 == 1024
    assert xsize1 == 512
    xmin = np.amin(x)
    xmax = np.amax(x)
    assert xmin == 4
    assert xmax == 515

    # test if channel 2 then right side of detector
    # test shape of x (1/2 detector) and min,max 516, 1027 (no reference pixels)
    input_model12 = datamodels.IFUImageModel((20, 20))
    input_model12.meta.instrument._instance.update(mirifushort_short)
    x, y = _get_2d_pixgrid(input_model12, '2')

    xsize0 = x.shape[0]
    xsize1 = x.shape[1]
    assert xsize0 == 1024
    assert xsize1 == 512
    xmin = np.amin(x)
    xmax = np.amax(x)
    assert xmin == 516
    assert xmax == 1027

    # test if channel 3 then right side of detector
    # test shape of x (1/2 detector) and min,max 516, 1027 (no reference pixels)
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    x, y = _get_2d_pixgrid(input_model34, '3')

    xsize0 = x.shape[0]
    xsize1 = x.shape[1]
    assert xsize0 == 1024
    assert xsize1 == 512
    xmin = np.amin(x)
    xmax = np.amax(x)
    assert xmin == 516
    assert xmax == 1027

    # test if channel 4 then left side of detector
    # test shape of x (1/2 detector) and min,max 4, 515 (includes reference pixels)
    input_model34 = datamodels.IFUImageModel((20, 20))
    input_model34.meta.instrument._instance.update(mirifulong_long)
    x, y = _get_2d_pixgrid(input_model34, '4')

    xsize0 = x.shape[0]
    xsize1 = x.shape[1]
    assert xsize0 == 1024
    assert xsize1 == 512
    xmin = np.amin(x)
    xmax = np.amax(x)
    assert xmin == 4
    assert xmax == 515
