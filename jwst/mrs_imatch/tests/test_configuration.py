"""
Unit test for Cube Build testing setting up configuration
"""

import pytest

from jwst import datamodels
from jwst.mrs_imatch.mrs_imatch_step import MRSIMatchStep
from jwst.mrs_imatch import mrs_imatch_step
import numpy as np

wcsinfo = {
    'dec_ref': -0.00244536159612126,
    'ra_ref': -0.00205553321270217,
    'roll_ref': -0.0,
    'v2_ref': -503.65447,
    'v3_ref': -318.74246,
    'v3yangle': 0.0,
    'vparity': -1
}

mirifushort_short = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'SHORT',
    'name': 'MIRI'
}

mirifushort_medium = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'MEDIUM',
    'name': 'MIRI'
}

mirifushort_long = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'LONG',
    'name': 'MIRI'
}

mirifulong_short = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'SHORT',
    'name': 'MIRI'
}

mirifulong_medium = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'MEDIUM',
    'name': 'MIRI'
}

mirifulong_long = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'LONG',
    'name': 'MIRI'
}


observation = {
    'date': '2019-01-01',
    'time': '17:00:00'}

subarray = {
    'fastaxis': 1,
    'name': 'FULL',
    'slowaxis': 2,
    'xsize': 1032,
    'xstart': 1,
    'ysize': 1024,
    'ystart': 1
}

@pytest.fixture(scope='function')
def miri_dither_ch12():
    """ Generate 4 dithered channel 12 data  """

    input_model1 = datamodels.IFUImageModel((20,20))
    input_model1.meta.wcsinfo._instance.update(wcsinfo)
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.observation._instance.update(observation)
    input_model1.meta.subarray._instance.update(subarray)
    input_model1.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model2 = datamodels.IFUImageModel()
    input_model2.meta.wcsinfo._instance.update(wcsinfo)
    input_model2.meta.instrument._instance.update(mirifushort_short)
    input_model2.meta.observation._instance.update(observation)
    input_model2.meta.subarray._instance.update(subarray)
    input_model2.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model3 = datamodels.IFUImageModel()
    input_model3.meta.wcsinfo._instance.update(wcsinfo)
    input_model3.meta.instrument._instance.update(mirifushort_short)
    input_model3.meta.observation._instance.update(observation)
    input_model3.meta.subarray._instance.update(subarray)
    input_model3.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model4 = datamodels.IFUImageModel()
    input_model4.meta.wcsinfo._instance.update(wcsinfo)
    input_model4.meta.instrument._instance.update(mirifushort_short)
    input_model4.meta.observation._instance.update(observation)
    input_model4.meta.subarray._instance.update(subarray)
    input_model4.meta.cal_step.assign_wcs = 'COMPLETE'

    # stuff in model container
    input_models = []
    input_models.append(input_model1)
    input_models.append(input_model2)
    input_models.append(input_model3)
    input_models.append(input_model4)

    return input_models


def test_imatch_background_test(_jail, miri_dither_ch12):
    """ Test if background subtract set or background exist correct thing is done"""

    all_models = datamodels.ModelContainer(miri_dither_ch12)
    
    # modify the data set backgroud subtracted
    new_container = []
    for m in all_models:
        m.meta.background.subtracted = True
        new_container.append(m)

    # test if background subtracted - raise error
    with pytest.raises(ValueError):
        step = MRSIMatchStep()
        step.run(new_container)


    # added a background and test is if it reset
    new_container = []
    degree = (1,1,1,)
    center = (5,5,5,)
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

    # test if background subtracted - raise error
    step = MRSIMatchStep()
    step._reset_background(new_container)
    
    for i in range(len(new_container)):
        m = new_container[i]
        test = len(m.meta.background.polynomial_info)
        assert test == 0
        

