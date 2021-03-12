import numpy as np
import numpy.random as rn

from jwst import datamodels
from jwst.barshadow import bar_shadow as bar


def test_create_shutter_elements():

    d1x1 = rn.random_sample((1001, 101))
    d1x3 = rn.random_sample((1001, 101))
    barshadow_model = datamodels.BarshadowModel(data1x1=d1x1, data1x3=d1x3)
    shutter_elements = bar.create_shutter_elements(barshadow_model)

    assert np.allclose(shutter_elements['first'], d1x1[:501, :], atol=1.e-10)
    assert np.allclose(shutter_elements['open_open'], d1x3[:501, :],
                       atol=1.e-10)
    assert np.allclose(shutter_elements['open_closed'], d1x1[500:, :],
                       atol=1.e-10)
    assert np.allclose(shutter_elements['closed_open'], d1x1[:501, :],
                       atol=1.e-10)
    assert np.allclose(shutter_elements['closed_closed'],
                       0.01 * np.ones((501, 101), dtype=np.float64),
                       atol=1.e-10)
    assert np.allclose(shutter_elements['last'], d1x1[501:, :], atol=1.e-10)


def test_create_first():

    shadow1x1 = rn.random_sample((1001, 101))

    shadow = bar.create_first(shadow1x1)

    assert np.allclose(shadow, shadow1x1[:501, :], atol=1.e-10)


def test_create_open_open():

    shadow1x3 = rn.random_sample((1001, 101))

    shadow = bar.create_open_open(shadow1x3)

    assert np.allclose(shadow, shadow1x3[:501, :], atol=1.e-10)


def test_create_open_closed():

    shadow1x1 = rn.random_sample((1001, 101))

    shadow = bar.create_open_closed(shadow1x1)

    assert np.allclose(shadow, shadow1x1[500:, :], atol=1.e-10)


def test_create_closed_open():

    shadow1x1 = rn.random_sample((1001, 101))

    shadow = bar.create_closed_open(shadow1x1)

    assert np.allclose(shadow, shadow1x1[:501, :], atol=1.e-10)


def test_create_closed_closed():

    shadow = bar.create_closed_closed()

    assert np.allclose(shadow, 0.01 * np.ones((501, 101), dtype=np.float64),
                       atol=1.e-10)


def test_create_last():

    shadow1x1 = rn.random_sample((1001, 101))

    shadow = bar.create_last(shadow1x1)

    assert np.allclose(shadow, shadow1x1[501:, :], atol=1.e-10)


def test_create_shadow():

    d1x1 = rn.random_sample((1001, 101))
    d1x3 = rn.random_sample((1001, 101))
    barshadow_model = datamodels.BarshadowModel(data1x1=d1x1, data1x3=d1x3)
    shutter_elements = bar.create_shutter_elements(barshadow_model)

    shutter_status = "11001x11"

    shadow = bar.create_shadow(shutter_elements, shutter_status)

    # first
    assert np.allclose(shadow[0:500, :], d1x1[0:500, :], atol=1.e-10)
    # open_open
    assert np.allclose(shadow[501:1000, :], d1x3[1:500, :], atol=1.e-10)
    # open_closed
    assert np.allclose(shadow[1001:1500, :], d1x1[501:1000, :], atol=1.e-10)
    # closed_closed
    assert np.allclose(shadow[1501:2000, :], 0.01, atol=1.e-10)
    # closed_open
    assert np.allclose(shadow[2001:2500, :], d1x1[1:500, :], atol=1.e-10)
    # open_open
    assert np.allclose(shadow[2501:3000, :], d1x3[1:500, :], atol=1.e-10)
    # open_open
    assert np.allclose(shadow[3001:3500, :], d1x3[1:500, :], atol=1.e-10)
    # open_open
    assert np.allclose(shadow[3501:4000, :], d1x3[1:500, :], atol=1.e-10)
    # last
    assert np.allclose(shadow[4001:4500, :], d1x1[502:1001, :], atol=1.e-10)


def test_create_empty_shadow_array():

    nshutters = 3
    shadow = bar.create_empty_shadow_array(nshutters)

    assert shadow.shape == (2000, 101)
    assert shadow.min() == 0.
    assert shadow.max() == 0.


def test_add_first_half_shutter():

    nshutters = 3
    shadow = bar.create_empty_shadow_array(nshutters)
    shadow_element = rn.random_sample((501, 101))
    shadow = bar.add_first_half_shutter(shadow, shadow_element)

    assert np.allclose(shadow[0:501, :], shadow_element, atol=1.e-10)


def test_add_next_shutter():

    nshutters = 3
    shadow = bar.create_empty_shadow_array(nshutters)
    shadow_element = rn.random_sample((501, 101))
    # The actual value of first_row would be a multiple of 500.
    first_row = 97
    dummy_value = 1.7

    shadow[:, :] = dummy_value
    shadow = bar.add_next_shutter(shadow, shadow_element, first_row)

    assert np.allclose(shadow[0:first_row, :], dummy_value, atol=1.e-10)

    # This column is the average of 1.7 and shadow_element[0, :].
    assert np.allclose(shadow[first_row, :],
                       (dummy_value + shadow_element[0, :]) / 2.,
                       atol=1.e-9)

    first_row += 1
    last_row = first_row + 500

    assert np.allclose(shadow[first_row:last_row, :], shadow_element[1:, :],
                       atol=1.e-10)

    assert np.allclose(shadow[last_row:, :], dummy_value, atol=1.e-10)


def test_add_last_half_shutter():

    nshutters = 7
    shadow = bar.create_empty_shadow_array(nshutters)
    shadow_element = rn.random_sample((501, 101))
    # The actual value of first_row would be a multiple of 500.
    first_row = 1819
    dummy_value = 1.3

    shadow[:, :] = dummy_value
    shadow = bar.add_last_half_shutter(shadow, shadow_element, first_row)

    assert np.allclose(shadow[0:first_row, :], dummy_value, atol=1.e-10)

    # This column is the average of 1.3 and shadow_element[0, :].
    assert np.allclose(shadow[first_row, :],
                       (dummy_value + shadow_element[0, :]) / 2.,
                       atol=1.e-9)

    first_row += 1
    last_row = first_row + 500

    assert np.allclose(shadow[first_row:last_row, :], shadow_element[1:, :],
                       atol=1.e-10)


def test_interpolate():

    d1x1 = np.arange(101 * 1001, dtype=np.float64) / (101. * 1001. - 1.)
    d1x1 = d1x1.reshape(1001, 101)
    d1x3 = d1x1.copy()
    barshadow_model = datamodels.BarshadowModel(data1x1=d1x1, data1x3=d1x3)
    shutter_elements = bar.create_shutter_elements(barshadow_model)

    shutter_status = "11x101"                   # 6 shutters
    # shadow will have shape (7 * 500, 101)     # 7 = len(shutter_status) + 1
    shadow = bar.create_shadow(shutter_elements, shutter_status)

    rows = np.arange(0, 3400 * 100 + 1, 10000, dtype=np.float64) / 100.
    columns = np.arange(0, 3400 * 100 + 1, 10000, dtype=np.float64) / 3400.
    rows = rows.reshape(5, 7)
    columns = columns.reshape(5, 7)

    correction = bar.interpolate(rows, columns, shadow, default=np.nan)

    compare = np.array([[0., 0.09993018, 0.19986036, 0.29979054,
                         0.39972072, 0.24989818, 0.10007564],
                        [0.20000582, 0.29993599, 0.39986617, 0.25004363,
                         0.10022110, 0.20015128, 0.30008147],
                        [0.40001165, 0.25018910, 0.10036656, 0.20029673,
                         0.30022691, 0.40015709, 0.50008730],
                        [0.60001744, 0.69994762, 0.79987779, 0.89980797,
                         0.50023272, 0.10065747, 0.20058765],
                        [0.30051783, 0.40044801, 0.50087771, 0.60130741,
                         0.70123758, 0.80116777, 0.90109789]])
    assert np.allclose(correction, compare, atol=1.e-6)


def test_has_uniform_source():

    data = np.zeros((10, 100), dtype=np.float32)
    slitlet = datamodels.SlitModel(data=data)

    # Since source_type has not been set yet, the barshadow step will
    # assume that the source is extended.
    assert bar.has_uniform_source(slitlet)

    slitlet.source_type = 'POINT'
    assert not bar.has_uniform_source(slitlet)          # not extended

    slitlet.source_type = 'UNKNOWN'
    # Since source_type is not 'POINT', the step will assume that the
    # source is extended.
    assert bar.has_uniform_source(slitlet)
