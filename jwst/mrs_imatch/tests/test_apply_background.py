"""
Unit test for mrs_imatch apply background
"""

import pytest

from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep
from jwst.mrs_imatch.mrs_imatch_step import apply_background_2d
import numpy as np

wcsinfo = {
    'dec_ref': 0.0,
    'ra_ref': 10.0,
    'roll_ref': 0.0,
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

    input_model1 = datamodels.IFUImageModel((30, 30))
    input_model1.meta.wcsinfo._instance.update(wcsinfo)
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.observation._instance.update(observation)
    input_model1.meta.subarray._instance.update(subarray)
    input_model1.meta.exposure.type = 'MIR_MRS'
    input_model1.data[:, :] = 10

    input_model2 = datamodels.IFUImageModel((30, 30))
    input_model2.meta.wcsinfo._instance.update(wcsinfo)
    input_model2.meta.instrument._instance.update(mirifushort_short)
    input_model2.meta.observation._instance.update(observation)
    input_model2.meta.subarray._instance.update(subarray)
    input_model2.meta.exposure.type = 'MIR_MRS'
    input_model2.data[:, :] = 10

    input_model3 = datamodels.IFUImageModel((30, 30))
    input_model3.meta.wcsinfo._instance.update(wcsinfo)
    input_model3.meta.instrument._instance.update(mirifushort_short)
    input_model3.meta.observation._instance.update(observation)
    input_model3.meta.subarray._instance.update(subarray)
    input_model3.meta.exposure.type = 'MIR_MRS'
    input_model3.data[:, :] = 20

    input_model4 = datamodels.IFUImageModel((30, 30))
    input_model4.meta.wcsinfo._instance.update(wcsinfo)
    input_model4.meta.instrument._instance.update(mirifushort_short)
    input_model4.meta.observation._instance.update(observation)
    input_model4.meta.subarray._instance.update(subarray)
    input_model4.meta.exposure.type = 'MIR_MRS'
    input_model4.data[:, :] = 10

    # stuff in model container
    input_models = []
    input_models.append(input_model1)
    input_models.append(input_model2)
    input_models.append(input_model3)
    input_models.append(input_model4)
    return input_models


def test_apply_background_2d(_jail, miri_dither_ch12):
    """ Test if  background polynomial is set it is subtracted correctly"""

    all_models = datamodels.ModelContainer(miri_dither_ch12)

    # test if given a background polynomial - apply_background_2d gives the correct answer
    # Polynomial information was created by running a faked full array data set to determine
    # what the polynomial values should be
    new_container = []
    degree = (1, 1, 1,)
    center = (9.99999408089876, -4.668612949274985e-06, 4.889999912818894)
    poly = np.ndarray(9)

    # run AssignWcsStep on first file and set wcs of other files to this

    channel = '1'
    for im, m in enumerate(all_models):
        if im == 0:
            m = AssignWcsStep.call(m)
            wcs1 = m.meta.wcs
            poly = [-2.50000000e+00, -6.66518331e-15,  5.67845589e-12, -1.23549218e-11,
                    3.26209108e-12, -5.43180357e-12, -2.54903452e-09,  6.21614553e-09]
        elif im == 1:
            m.meta.wcs = wcs1
            poly = [-2.50000000e+00, -9.29031986e-15,  1.74437380e-12, -3.94894956e-12,
                     4.71729481e-13,  5.42031845e-13, -9.96151554e-10,  2.78281950e-09]
        elif im == 2:
            m.meta.wcs = wcs1
            poly = [7.50000000e+00,  2.22921836e-14, -6.97131279e-12,  1.58336906e-11,
                     -5.43704212e-12,  6.79109678e-12,  2.63515372e-09, -8.00226976e-09]
        else:
            m.meta.wcs = wcs1
            poly = [-2.50000000e+00, -6.33668043e-15, -4.51516905e-13,  4.70180795e-13,
                     1.70322157e-12, -1.90132506e-12,  9.10032350e-10, -9.96695272e-10]

        poly = np.asarray(poly)
        m.meta.background.polynomial_info.append(
            {
                'degree': degree,
                'refpoint': center,
                'coefficients': poly.ravel().tolist(),
                'channel': channel
                }
            )

        # input 10 output 12.5 for all images
        apply_background_2d(m, channel, subtract=True)
        new_container.append(m)

    # all the data should now be the same for the science pixels
    # the edge pixels are located in non-science region and no
    # background subtraction is done on these pixels

    data1 = new_container[0].data[:, 16:]
    data2 = new_container[1].data[:, 16:]
    data3 = new_container[2].data[:, 16:]
    data4 = new_container[3].data[:, 16:]

    assert(np.allclose(data1, data2, rtol=1e-6))
    assert(np.allclose(data2, data3, rtol=1e-6))
    assert(np.allclose(data3, data4, rtol=1e-6))
    assert(np.allclose(data1, data4, rtol=1e-6))
