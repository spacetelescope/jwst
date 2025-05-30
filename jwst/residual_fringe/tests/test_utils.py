import numpy as np
import pytest

from jwst.residual_fringe import utils


@pytest.fixture()
def slice_map():
    """
    Make a slice map image containing 4 slices.

    There will be two channel 1 slices and two channel 2.
    """
    shape = (20, 90)
    map_image = np.zeros(shape)
    map_image[:, 10:20] = 101
    map_image[:, 30:40] = 102
    map_image[:, 50:60] = 201
    map_image[:, 70:80] = 202
    return map_image


def test_slice_info_ch1(slice_map):
    result = utils.slice_info(slice_map, 1)
    slices_in_channel, xrange_channel, slice_x_ranges, all_slice_masks = result
    assert np.all(slices_in_channel == [101, 102])
    assert np.all(xrange_channel == [9, 40])
    assert np.all(slice_x_ranges == [[101, 9, 20], [102, 29, 40]])
    assert np.sum(all_slice_masks[0]) == 200
    assert np.sum(all_slice_masks[1]) == 200


def test_slice_info_ch2(slice_map):
    result = utils.slice_info(slice_map, 2)
    slices_in_channel, xrange_channel, slice_x_ranges, all_slice_masks = result
    assert np.all(slices_in_channel == [201, 202])
    assert np.all(xrange_channel == [49, 80])
    assert np.all(slice_x_ranges == [[201, 49, 60], [202, 69, 80]])
    assert np.sum(all_slice_masks[0]) == 200
    assert np.sum(all_slice_masks[1]) == 200
