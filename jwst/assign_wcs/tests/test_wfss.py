import os.path
import numpy as np
from numpy.testing import assert_allclose

from astropy.modeling.models import (Polynomial1D, Polynomial2D, Shift,
                                     Const1D, Mapping)
from gwcs import wcs
from gwcs.wcstools import grid_from_bounding_box

from stdatamodels.jwst.datamodels import SlitModel
from stdatamodels.jwst.transforms import models as transforms

from jwst.extract_2d.grisms import compute_wavelength_array

import pytest

from jwst.assign_wcs.tests import data
data_path = os.path.split(os.path.abspath(data.__file__))[0]


def create_nircam_slit(model, x0, y0, order):
    """ Create a SlitModel representing a grism slit."""
    ymin = 0
    xmin = 0
    # ymax = 58
    # xmax = 1323
    model = Mapping((0, 1, 0, 0, 0)) | (Shift(xmin) & Shift(ymin) &
                                        Const1D(x0) & Const1D(y0) & Const1D(order)) | model
    wcsobj = wcs.WCS([('det', model), ('world', None)])
    wcsobj.bounding_box = ((20, 25), (800, 805))
    slit = SlitModel()
    slit.meta.wcs = wcsobj
    slit.source_xpos = x0
    slit.source_ypos = y0
    return slit


def create_niriss_slit(model, x0, y0, order):
    """ Create a SlitModel representing a grism slit."""
    ymin = y0
    xmin = x0
    model = Mapping((0, 1, 0, 0, 0)) | (Shift(xmin) & Shift(ymin) &
                                        Const1D(x0) & Const1D(y0) &
                                        Const1D(order)) | model
    wcsobj = wcs.WCS([('det', model), ('world', None)])
    slit = SlitModel()
    slit.meta.wcs = wcsobj
    slit.source_xpos = x0
    slit.source_ypos = y0
    return slit


@pytest.fixture(scope='module')
def niriss_models():
    ymodels = [[Polynomial2D(2, c0_0=59.70176, c1_0=-0.00000036, c2_0=0., c0_1=-0.00000026),
                Polynomial2D(2, c0_0=-330.4911, c1_0=0.00000089, c2_0=-0., c0_1=0.00000064),
                Polynomial2D(2, c0_0=0.00041831, c1_0=-0.00000055, c2_0=0., c0_1=-0.0000004)],
               [Polynomial2D(2, c0_0=-100.4272, c1_0=-0.00000005, c2_0=0., c0_1=-0.0000004),
                Polynomial2D(2, c0_0=-662.3936, c1_0=0.00000005, c2_0=-0., c0_1=0.00000103),
                Polynomial2D(2, c0_0=0.00025086, c1_0=0.00000003, c2_0=0., c0_1=-0.00000066)],
               [Polynomial2D(2, c0_0=-260.8984, c1_0=0.00000052, c2_0=0., c0_1=0.0000013),
                Polynomial2D(2, c0_0=-993.5874, c1_0=-0.00000134, c2_0=-0., c0_1=-0.00000336),
                Polynomial2D(2, c0_0=-0.00152898, c1_0=0.00000086, c2_0=0., c0_1=0.00000215)],
               [Polynomial2D(2, c0_0=379.5311, c1_0=-0.00000027, c2_0=0., c0_1=-0.00000025),
                Polynomial2D(2, c0_0=330.49, c1_0=0.00000072, c2_0=-0., c0_1=0.00000066),
                Polynomial2D(2, c0_0=0.00025935, c1_0=-0.00000048, c2_0=0., c0_1=-0.00000044)],
               [Polynomial2D(2, c0_0=217.1429, c1_0=-0.00000001),
                Polynomial2D(2, c0_0=-7.38096, c1_0=0.00000002, c0_1=-0.00000001),
                Polynomial2D(2, c0_0=0.00000253, c1_0=-0.00000001)]]

    xmodels = [[Polynomial2D(2, c0_0=-0.2846298),
                Polynomial2D(2, c0_0=-0.08449617),
                Polynomial2D(2, c0_0=-0.00243989)],
               [Polynomial2D(2, c0_0=-0.412634),
                Polynomial2D(2, c0_0=-0.08853529),
                Polynomial2D(2, c0_0=0.0004952)],
               [Polynomial2D(2, c0_0=0.09586888),
                Polynomial2D(2, c0_0=-0.5574736),
                Polynomial2D(2, c0_0=-0.00080743)],
               [Polynomial2D(2, c0_0=0.3766071),
                Polynomial2D(2, c0_0=0.3271358),
                Polynomial2D(2, c0_0=-0.00487867)],
               [Polynomial2D(2, c0_0=-0.99),
                Polynomial2D(2),
                Polynomial2D(2)]]

    lmodels = [Polynomial1D(1, c0=0.75, c1=1.55),
               Polynomial1D(1, c0=0.75, c1=1.55),
               Polynomial1D(1, c0=0.75, c1=1.55),
               Polynomial1D(1, c0=0.75, c1=1.55),
               Polynomial1D(1, c0=0.75, c1=1.55)]
    return xmodels, ymodels, lmodels


def test_NIRCAMForwardRowGrismDispersion():
    xmodels = [Polynomial1D(1, c0=0.59115385, c1=0.00038615),
               Polynomial1D(1, c0=-0.16596154, c1=0.00019308)]
    ymodels = [Polynomial1D(1, c0=0., c1=0.), Polynomial1D(1, c0=0., c1=0.)]
    lmodels = [Polynomial1D(1, c0=2.4, c1=2.6), Polynomial1D(1, c0=2.4, c1=2.6)]
    model = transforms.NIRCAMForwardRowGrismDispersion([1, 2], lmodels, xmodels, ymodels)

    x0 = 913.7
    y0 = 15.5
    order = 1

    slit = create_nircam_slit(model, x0, y0, order)

    expected = np.array([[3.03973415, 3.04073814, 3.04174213, 3.04274612, 3.04375011, 3.0447541],
                         [3.03973415, 3.04073814, 3.04174213, 3.04274612, 3.04375011, 3.0447541],
                         [3.03973415, 3.04073814, 3.04174213, 3.04274612, 3.04375011, 3.0447541],
                         [3.03973415, 3.04073814, 3.04174213, 3.04274612, 3.04375011, 3.0447541],
                         [3.03973415, 3.04073814, 3.04174213, 3.04274612, 3.04375011, 3.0447541],
                         [3.03973415, 3.04073814, 3.04174213, 3.04274612, 3.04375011, 3.0447541]])

    # refactored call
    x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
    wavelength = compute_wavelength_array(slit)  # x, y, np.zeros(x.shape)  +x0, np.zeros(y.shape)+y0, np.zeros(x.shape)+order)
    assert_allclose(wavelength, expected)

    with pytest.raises(ValueError):
        slit = create_nircam_slit(model, x0, y0, 3)
        compute_wavelength_array(slit)


def test_NIRCAMForwardColumnGrismDispersion():
    ymodels = [Polynomial1D(1, c0=0.5911538461431823, c1=0.000386153846153726),
               Polynomial1D(1, c0=-0.1659615384582264, c1=0.0001930769230768787)]

    xmodels = [Polynomial1D(1, c0=0., c1=0.), Polynomial1D(1, c0=0., c1=0.)]
    lmodels = [Polynomial1D(1, c0=2.4, c1=2.6), Polynomial1D(1, c0=2.4, c1=2.6)]
    model = transforms.NIRCAMForwardColumnGrismDispersion([1, 2], lmodels, xmodels, ymodels)

    x0 = 913.7
    y0 = 15.5
    order = 1

    slit = create_nircam_slit(model, x0, y0, order)
    expected = np.array([[4.724638, 4.724638, 4.724638, 4.724638, 4.724638, 4.724638],
                         [4.725642, 4.725642, 4.725642, 4.725642, 4.725642, 4.725642],
                         [4.726646, 4.726646, 4.726646, 4.726646, 4.726646, 4.726646],
                         [4.72765, 4.72765, 4.72765, 4.72765, 4.72765, 4.72765],
                         [4.728654, 4.728654, 4.728654, 4.728654, 4.728654, 4.728654],
                         [4.729658, 4.729658, 4.729658, 4.729658, 4.729658, 4.729658]])

    # refactored call
    x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
    wavelength = compute_wavelength_array(slit)  # x, y, np.zeros(x.shape)  +x0, np.zeros(y.shape)+y0, np.zeros(x.shape)+order)
    assert_allclose(wavelength, expected)

    with pytest.raises(ValueError):
        slit = create_nircam_slit(model, x0, y0, 3)
        compute_wavelength_array(slit)


def test_NIRCAMBackwardDispersion():
    forward_ymodels = [Polynomial1D(1, c0=0.5911538461431823, c1=0.000386153846153726),
                       Polynomial1D(1, c0=-0.1659615384582264, c1=0.0001930769230768787)]

    forward_xmodels = [Polynomial1D(1, c0=0., c1=0.), Polynomial1D(1, c0=0., c1=0.)]
    forward_lmodels = [Polynomial1D(1, c0=2.4, c1=2.6), Polynomial1D(1, c0=2.4, c1=2.6)]

    forward_model = transforms.NIRCAMForwardColumnGrismDispersion([1, 2], forward_lmodels,
                                                                  forward_xmodels, forward_ymodels)

    ymodels = [Polynomial1D(1, c0=-1530.8764939967652, c1=2589.641434263754),
               Polynomial1D(1, c0=859.5617529710912, c1=5179.282868527087)]

    xmodels = [Polynomial1D(1, c0=0., c1=0.),
               Polynomial1D(1, c0=0., c1=0.)]
    lmodels = [Polynomial1D(1, c0=-0.923076923076923, c1=0.3846153846153846),
               Polynomial1D(1, c0=-0.923076923076923, c1=0.3846153846153846)]
    model = transforms.NIRCAMBackwardGrismDispersion([1, 2], lmodels, xmodels, ymodels)
    wavelength = np.array([[4.724638, 4.724638, 4.724638, 4.724638, 4.724638, 4.724638],
                           [4.725642, 4.725642, 4.725642, 4.725642, 4.725642, 4.725642],
                           [4.726646, 4.726646, 4.726646, 4.726646, 4.726646, 4.726646],
                           [4.72765, 4.72765, 4.72765, 4.72765, 4.72765, 4.72765],
                           [4.728654, 4.728654, 4.728654, 4.728654, 4.728654, 4.728654],
                           [4.729658, 4.729658, 4.729658, 4.729658, 4.729658, 4.729658]])

    x0 = 913.7
    y0 = 15.5
    order = 1

    slit = create_nircam_slit(forward_model, x0, y0, order)

    expected_xdx = np.array([[20., 21., 22., 23., 24., 25.],
                             [20., 21., 22., 23., 24., 25.],
                             [20., 21., 22., 23., 24., 25.],
                             [20., 21., 22., 23., 24., 25.],
                             [20., 21., 22., 23., 24., 25.],
                             [20., 21., 22., 23., 24., 25.]])

    expected_ydy = np.array([[1584.50000003, 1584.50000003, 1584.50000003, 1584.50000003, 1584.50000003, 1584.50000003],
                             [1586.50000003, 1586.50000003, 1586.50000003, 1586.50000003, 1586.50000003, 1586.50000003],
                             [1588.50000003, 1588.50000003, 1588.50000003, 1588.50000003, 1588.50000003, 1588.50000003],
                             [1590.50000003, 1590.50000003, 1590.50000003, 1590.50000003, 1590.50000003, 1590.50000003],
                             [1592.50000003, 1592.50000003, 1592.50000003, 1592.50000003, 1592.50000003, 1592.50000003],
                             [1594.50000003, 1594.50000003, 1594.50000003, 1594.50000003, 1594.50000003, 1594.50000003]])
    # refactored call
    x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
    xdx, ydy, _, _, _ = model(x, y, wavelength, np.zeros(x.shape) + order)
    assert_allclose(xdx, expected_xdx)
    assert_allclose(ydy, expected_ydy)


def test_NIRISSBackwardDispersion(niriss_models):
    forward_xmodels, forward_ymodels, forward_lmodels = niriss_models

    forward_model = transforms.NIRISSForwardRowGrismDispersion([1, 2, 3, -1], forward_lmodels,
                                                               forward_xmodels, forward_ymodels,
                                                               theta=0)

    # NirissBackward model uses xmodels, ymodels and invlmodels
    lmodels = [Polynomial1D(1, c0=-0.48387097, c1=0.64516129),
               Polynomial1D(1, c0=-0.48387097, c1=0.64516129),
               Polynomial1D(1, c0=-0.48387097, c1=0.64516129),
               Polynomial1D(1, c0=-0.48387097, c1=0.64516129),
               Polynomial1D(1, c0=-0.48387097, c1=0.64516129)]

    model = transforms.NIRISSBackwardGrismDispersion([1, 2, 3, -1], lmodels=lmodels,
                                                     xmodels=forward_xmodels,
                                                     ymodels=forward_ymodels,
                                                     theta=0.0)

    wavelength = np.array(
        [[2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3],
         [0.98553179, 0.98553179, 0.98553179, 0.98553179, 0.98553179, 0.98553179, 0.98553179],
         [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75],
         [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75],
         [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75],
         [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75],
         [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75]])

    x0 = 913.7
    y0 = 15.5
    order = 1

    slit = create_niriss_slit(forward_model, x0, y0, order)
    slit.meta.wcs.bounding_box = ((910, 916), (12, 18))

    expected_xdx = np.array(
        [[909.62843414, 910.62843414, 911.62843414, 912.62843414,
          913.62843414, 914.62843414, 915.62843414],
         [909.70247416, 910.70247416, 911.70247416, 912.70247416,
          913.70247416, 914.70247416, 915.70247416],
         [909.7153702, 910.7153702, 911.7153702, 912.7153702,
          913.7153702, 914.7153702, 915.7153702],
         [909.7153702, 910.7153702, 911.7153702, 912.7153702,
          913.7153702, 914.7153702, 915.7153702],
         [909.7153702, 910.7153702, 911.7153702, 912.7153702,
          913.7153702, 914.7153702, 915.7153702],
         [909.7153702, 910.7153702, 911.7153702, 912.7153702,
          913.7153702, 914.7153702, 915.7153702],
         [909.7153702, 910.7153702, 911.7153702, 912.7153702,
          913.7153702, 914.7153702, 915.7153702]])

    expected_ydy = np.array(
        [[-258.78893914, -258.78893916, -258.78893918, -258.7889392,
          -258.78893922, -258.78893924, -258.78893926],
         [22.48144873, 22.48144849, 22.48144825, 22.48144802,
          22.48144778, 22.48144754, 22.4814473],
         [73.70142959, 73.70142923, 73.70142887, 73.70142851,
          73.70142815, 73.70142779, 73.70142743],
         [74.70142933, 74.70142897, 74.70142861, 74.70142825,
          74.70142789, 74.70142753, 74.70142717],
         [75.70142907, 75.70142871, 75.70142835, 75.70142799,
          75.70142763, 75.70142727, 75.70142691],
         [76.70142881, 76.70142845, 76.70142809, 76.70142773,
          76.70142737, 76.70142701, 76.70142665],
         [77.70142855, 77.70142819, 77.70142783, 77.70142747,
          77.70142711, 77.70142675, 77.70142639]])

    # refactored call
    x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
    xdx, ydy, _, _, _ = model(x, y, wavelength, np.zeros(x.shape) + 1)
    assert_allclose(xdx, expected_xdx)
    assert_allclose(ydy, expected_ydy)


def test_NIRISSForwardRowGrismDispersion(niriss_models):
    xmodels, ymodels, lmodels = niriss_models
    model = transforms.NIRISSForwardRowGrismDispersion([1, 2, 3, -1], lmodels, xmodels,
                                                       ymodels, theta=354.222)

    x0 = 913.7
    y0 = 15.5
    order = 1

    slit = create_niriss_slit(model, x0, y0, order)
    slit.meta.wcs.bounding_box = ((910, 916), (12, 18))

    wavelength = compute_wavelength_array(slit)
    expected = np.array(
        [[-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588],
         [-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588],
         [-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588],
         [-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588],
         [-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588],
         [-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588],
         [-41.26984722, -41.31631533, -41.36278344, -41.40925155,
          -41.45571966, -41.50218777, -41.54865588]])
    # refactored call
    x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
    wavelength = compute_wavelength_array(slit)
    assert_allclose(wavelength, expected)

#
# def test_NIRISSForwardColumnGrismDispersion(niriss_models):
#     xmodels, ymodels, lmodels = niriss_models
#
#     model = transforms.NIRISSForwardColumnGrismDispersion([1, 2, 3, -1], lmodels=lmodels,
#                                                           xmodels=xmodels, ymodels=ymodels,
#                                                           theta=33.5677)
#
#     x0 = 913.7
#     y0 = 15.5
#     order = 1
#
#     slit = create_niriss_slit(model, x0, y0, order)
#     slit.meta.wcs.bounding_box = ((910, 916), (12, 18))
#
#     expected = np.array(
#         [[1.05844596, 1.05844596, 1.05844596, 1.05844596, 1.05844596,
#           1.05844596, 1.05844596],
#          [1.0500404, 1.0500404, 1.0500404, 1.0500404, 1.0500404,
#             1.0500404, 1.0500404],
#             [1.04163483, 1.04163483, 1.04163483, 1.04163483, 1.04163483,
#              1.04163483, 1.04163483],
#             [1.03322927, 1.03322927, 1.03322927, 1.03322927, 1.03322927,
#              1.03322927, 1.03322927],
#             [1.02482371, 1.02482371, 1.02482371, 1.02482371, 1.02482371,
#              1.02482371, 1.02482371],
#             [1.01641815, 1.01641815, 1.01641815, 1.01641815, 1.01641815,
#              1.01641815, 1.01641815],
#             [1.00801258, 1.00801258, 1.00801258, 1.00801258, 1.00801258,
#              1.00801258, 1.00801258]])
#
#     # refactored call
#     x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
#     wavelength = compute_wavelength_array(slit)
#     assert_allclose(wavelength, expected)
