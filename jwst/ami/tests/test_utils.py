"""Unit tests for AMI utils module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
import stdatamodels.jwst.datamodels as dm

from jwst.ami import utils
from jwst.ami.bp_fix import filthp_d, filtwl_d


def test_affine2d():
    pass


@pytest.mark.parametrize("n", [20, 21])
def test_makedisk(n):
    r = 4
    arr = utils.makedisk(n, r)

    # test disk is symmetric about center in both directions
    assert_allclose(arr[:10, :], arr[-10:, :][::-1, :])
    assert_allclose(arr[:, :10], arr[:, -10:][:, ::-1])

    # test disk is all 0 or 1
    assert np.all(np.logical_or(arr == 0, arr == 1))


@pytest.mark.parametrize("dtype", [int, float])
def test_avoidhexsingularity(dtype):
    rotation = np.arange(-91, 92, 1, dtype=dtype)
    for rot in rotation:
        is_multiple_of_15 = rot % 15 == 0
        rot_fixed = utils.avoidhexsingularity(rot)
        if is_multiple_of_15:
            assert rot_fixed != rot
        else:
            assert rot_fixed == rot


@pytest.mark.parametrize(
    "shape, center",
    [
        ((10, 10), (4.5, 4.5)),
        ((11, 11), (5, 5)),
    ],
)
def test_centerpoint(shape, center):
    assert utils.centerpoint(shape) == center


def test_min_distance_to_edge():
    """
    cntrimg=True case still needs testing.
    But this will not be used if center_imgpeak can be deleted.
    """
    sz = 21
    shp = (sz, sz)
    img = np.zeros(shp)
    img[10, 10] = 1

    # basic test
    x_out, y_out, dist = utils.min_distance_to_edge(img)
    assert dist == 10.0
    assert x_out == 10
    assert y_out == 10

    # add an equally bright source closer to the edge
    img[5, 8] = 1
    x_out, y_out, dist = utils.min_distance_to_edge(img)
    assert dist == 5.0
    assert x_out == 5
    assert y_out == 8

    # add a fainter source closer to the edge
    img[4, 8] = 0.5
    x_out, y_out, dist = utils.min_distance_to_edge(img)
    assert dist == 5.0
    assert x_out == 5
    assert y_out == 8

    # No failure when no sources, just dist goes to zero
    img = np.zeros(shp)
    x_out, y_out, dist = utils.min_distance_to_edge(img)
    assert dist == 0.0


def test_find_centroid():
    """Tests of find_centroid and findslope functions."""

    arr = np.zeros((30, 30), dtype="f4")
    arr[17, 15] = 1
    htilt, vtilt = utils.find_centroid(arr)
    assert_allclose(htilt, 2.5)
    assert_allclose(vtilt, 0.5)


def test_quadratic_extremum():
    """Tests of quadratic_extremum and findpeak_1d functions"""

    def quadratic(p, x):
        """Ordering of p is consistent with np.polyfit"""
        return p[0] * x**2 + p[1] * x + p[2]

    x = np.linspace(-5, 5, 100)

    # Test of quadratic_extremum with a maximum
    p = np.array([-2.0, 3.0, 7.0])
    maxx, maxy = utils.quadratic_extremum(p)
    true_maxx = 0.75
    true_maxy = 8.125
    assert_allclose([maxx, maxy], [true_maxx, true_maxy])

    # test findpeak_1d will find the same maximum
    y = quadratic(p, x)
    maxx, maxy = utils.findpeak_1d(x, y)
    assert_allclose([maxx, maxy], [true_maxx, true_maxy])

    # Test of quadratic_extremum with a minimum
    p = np.array([2.0, 3.0, 7.0])
    minx, miny = utils.quadratic_extremum(p)
    true_minx = -0.75
    true_miny = 5.875
    assert_allclose([minx, miny], [true_minx, true_miny])

    # test findpeak_1d will find the same minimum
    y = quadratic(p, x)
    minx, miny = utils.findpeak_1d(x, y)
    assert_allclose([minx, miny], [true_minx, true_miny])


def test_make_a():
    """Test of make_a in utils module"""
    nh = 4  # number of holes
    arr = utils.make_a(nh)

    true_arr = np.array(
        [
            [-1.0, 1.0, 0.0, 0.0],
            [-1.0, 0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0, 1.0],
            [0.0, -1.0, 1.0, 0.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, -1.0, 1.0],
        ]
    )
    assert_allclose(arr, true_arr)


def test_fringes2pistons():
    """Test of fringes2pistons in utils module"""
    fringephases = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    nholes = 5

    result = utils.fringes2pistons(fringephases, nholes)

    true_result = np.array([-0.02, -0.034, -0.02, 0.014, 0.06])
    assert_allclose(result, true_result)


def test_rebin():
    """Test of rebin() and krebin() in utils module"""
    arr = np.arange(24).reshape((3, 8)) / 10.0
    rc = tuple((2, 2))

    binned_arr = utils.rebin(arr, rc)

    true_arr = np.array([[5.1, 6.3, 7.5, 8.7]])
    assert_allclose(binned_arr, true_arr)


def test_rcrosscorrelate():
    """Test of rcrosscorrelate() in utils module"""
    a = np.array(
        [
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0, 13.0],
            [14.0, 15.0, 16.0, 17.0],
        ]
    )

    b = np.array(
        [
            [-5.0, -4.0, -3.0, -2.0],
            [-1.0, 0.0, 1.0, 2.0],
            [3.0, 4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0, 10.0],
        ]
    )

    result = utils.rcrosscorrelate(a, b)

    true_result = np.array(
        [
            [0.19865015, 0.20767971, 0.23476836, 0.20767971],
            [0.34312299, 0.35215254, 0.3792412, 0.35215254],
            [0.77654151, 0.78557106, 0.81265972, 0.78557106],
            [0.34312299, 0.35215254, 0.3792412, 0.35215254],
        ]
    )
    assert_allclose(result, true_result)


def test_crosscorrelate():
    """Test of crosscorrelate() in utils module"""
    a = np.array(
        [
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0, 13.0],
            [14.0, 15.0, 16.0, 17.0],
        ]
    )

    b = np.array(
        [
            [-5.0, -4.0, -3.0, -2.0],
            [-1.0, 0.0, 1.0, 2.0],
            [3.0, 4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0, 10.0],
        ]
    )

    result = utils.crosscorrelate(a, b)

    true_result = np.array(
        [
            [176.0 + 0.0j, 184.0 + 0.0j, 208.0 + 0.0j, 184.0 + 0.0j],
            [304.0 + 0.0j, 312.0 + 0.0j, 336.0 + 0.0j, 312.0 + 0.0j],
            [688.0 + 0.0j, 696.0 + 0.0j, 720.0 + 0.0j, 696.0 + 0.0j],
            [304.0 + 0.0j, 312.0 + 0.0j, 336.0 + 0.0j, 312.0 + 0.0j],
        ]
    )
    assert_allclose(result, true_result)


def test_rotate2dccw():
    # test 90 degree rotation
    vectors = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    rotated = utils.rotate2dccw(vectors, np.pi / 2)
    # very near zero so default rtol
    assert_allclose(rotated, np.array([[0, 1], [-1, 0], [0, -1], [1, 0]]), atol=1e-10)

    # test that small positive rotation makes x smaller and y larger
    vector = np.array(
        [
            [1, 1],
        ]
    )
    rotated = utils.rotate2dccw(vector, np.pi / 30)[0]
    assert rotated[0] < vector[0][0]
    assert rotated[1] > vector[0][1]


def test_get_cw_beta(example_model, bandpass):
    filt = example_model.meta.instrument.filter
    cw, beta = utils.get_cw_beta(bandpass)

    # we constructed the bandpass to be symmetric around the input filter central wave
    # so the central wavelength we get here should equal that input
    assert np.isclose(cw, filtwl_d[filt])

    # similarly, we should get back the same fractional width as we would get from
    # the input half-power points, which are where the top-hat function cuts
    # on and cuts off. This is affected by the wavelength sampling so the rtol
    # being relatively larger (1e-3) is not a concern.
    wl_min, wl_max = filthp_d[filt]
    assert np.isclose(beta, (wl_max - wl_min) / cw, rtol=1e-3)


def test_handle_bandpass(bandpass, throughput_model):
    """Tests for get_filt_spec, get_flat_spec, combine_src_filt, and handle_bandpass functions"""

    # test case where user-defined bandpass is None and we are relying on the throughput model
    bandpass0 = utils.handle_bandpass(None, throughput_model)
    wls = bandpass0[:, 1]
    thrus = bandpass0[:, 0]

    # check data type
    for arr in [wls, thrus]:
        assert isinstance(arr, np.ndarray)
        # this gets upcast to float64 for some reason, even though it's specified to be
        # float32 by the throughput_model
        # assert arr.dtype == np.dtype("float32")
        assert arr.size > 1

    # ensure the wavelengths are in meters and all within the half-power points
    # cutoff to half-power points is enforced by the fact that this is where throughput model
    # switches from equaling unity to equaling 1e-3
    filt = throughput_model.meta.instrument.filter
    wlhp_l, wlhp_h = filthp_d[filt]
    assert np.all(wls >= wlhp_l)
    assert np.all(wls <= wlhp_h)

    # ensure the throughput is flat and normalized
    assert np.allclose(thrus / np.mean(thrus), 1.0, rtol=1e-5)
    assert np.isclose(np.sum(thrus), 1.0)

    # test case where user-defined bandpass is provided as an array
    # In this case, the bandpass just gets returned without modification
    bandpass1 = utils.handle_bandpass(bandpass, throughput_model)
    assert np.all(bandpass1 == bandpass)

    # test case where user-defined bandpass is provided as a synphot object
    # Again the bandpass just gets converted to array then returned without modification
    bandpass_synphot = utils.get_filt_spec(throughput_model)
    bandpass2 = utils.handle_bandpass(bandpass_synphot, throughput_model)

    # use atol here because bandpass is normalized to unity.
    # There is some small numerical problem that makes some of the low-throughput values
    # not exactly equal to 0.001 in the synphot object, but the atol is small enough
    # that we don't care
    assert np.allclose(bandpass2, bandpass1, atol=1e-4)


def test_degrees_per_pixel(example_model):
    # example_model has no wcsinfo to start with
    pixscale_x, pixscale_y = utils.degrees_per_pixel(example_model)
    assert pixscale_x == pixscale_y == 65.6 / (60.0 * 60.0 * 1000)

    # add cdelt to example_model
    example_model.meta.wcsinfo.cdelt1 = 0.1
    example_model.meta.wcsinfo.cdelt2 = 0.05
    pixscale_x, pixscale_y = utils.degrees_per_pixel(example_model)
    assert pixscale_x == 0.1
    assert pixscale_y == 0.05

    # add CD matrix, which should take precedence over cdelt
    cd = np.array([[0.2, 0.0], [0.0, 0.1]])
    example_model.meta.wcsinfo.cd1_1 = cd[0, 0]
    example_model.meta.wcsinfo.cd1_2 = cd[0, 1]
    example_model.meta.wcsinfo.cd2_1 = cd[1, 0]
    example_model.meta.wcsinfo.cd2_2 = cd[1, 1]
    pixscale_x, pixscale_y = utils.degrees_per_pixel(example_model)
    assert pixscale_x == 0.2
    assert pixscale_y == 0.1
