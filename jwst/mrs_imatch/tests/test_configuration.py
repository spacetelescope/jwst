"""
Unit test for Cube Build testing setting up configuration
"""

import pytest

from jwst import datamodels
from jwst.mrs_imatch.mrs_imatch_step import MRSIMatchStep
import numpy as np

mirifushort_short = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'SHORT',
    'name': 'MIRI'
}


@pytest.fixture(scope='function')
def miri_dither_ch12():
    """ Generate 4 dithered channel 12 data  """

    input_model1 = datamodels.IFUImageModel((20,20))
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model2 = datamodels.IFUImageModel()
    input_model2.meta.instrument._instance.update(mirifushort_short)
    input_model2.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model3 = datamodels.IFUImageModel()
    input_model3.meta.instrument._instance.update(mirifushort_short)
    input_model3.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model4 = datamodels.IFUImageModel()
    input_model4.meta.instrument._instance.update(mirifushort_short)
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


def test_imatch_degree_test(_jail, miri_dither_ch12):
    """ Test if polynomial degree is not a 3 element tuple or integer - raise error """

    all_models = datamodels.ModelContainer(miri_dither_ch12)
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

    # test catches degree set incorrectly - raise error
    #  check that degree must be a 3 element tuple
    with pytest.raises(ValueError):
        step = MRSIMatchStep()
        step.bkg_degree = (1,1,)
        step.run(new_container)
