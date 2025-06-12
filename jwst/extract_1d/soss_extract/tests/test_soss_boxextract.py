import pytest
import numpy as np

from .conftest import DATA_SHAPE
from jwst.extract_1d.soss_extract.soss_boxextract import (
    get_box_weights,
    box_extract,
    estim_error_nearest_data,
)


WIDTH = 5.1


@pytest.fixture()
def box_weights(trace1d):
    weights_list = []
    for order in [0, 1]:
        tracex, tracey, wavetrace = trace1d[order]
        weights_list.append(get_box_weights(tracey, WIDTH, DATA_SHAPE, tracex))
    return weights_list


def test_get_box_weights(trace1d, box_weights):
    """
    Order 1 is easy because tracex.size and tracey.size are equal to data_shape[1]
    Order 2 tests the case where they are not equal
    """

    for order in [0, 1]:
        tracex, tracey, wavetrace = trace1d[order]
        weights = box_weights[order]

        # check weights are between zero and 1
        assert weights.shape == DATA_SHAPE
        assert np.max(weights) == 1
        assert np.min(weights) == 0

        weight_sum = np.sum(weights, axis=0)
        # some weights are zero because of the trace profile cutoff in order 2
        assert np.sum(weight_sum == 0) == (DATA_SHAPE[1] - tracey.size)

        # check sum of weights across y-axis is width
        weight_sum = weight_sum[np.nonzero(weight_sum)]
        assert np.allclose(weight_sum, WIDTH)

        # check at least some partial weights exist
        assert np.sum(weight_sum == 1) < weight_sum.size


def test_box_extract(trace1d, box_weights, imagemodel):
    data, err = imagemodel
    mask = np.isnan(data)

    for order in [0, 1]:
        weights = box_weights[order]
        cols, flux, flux_err, npix = box_extract(data, err, mask, weights)

        # test that cols just represents the data
        assert np.allclose(cols, np.arange(data.shape[1]))

        # test flux and flux_err are NaN where order 2 is cut off, but have good values elsewhere
        xtrace = trace1d[order][0]
        for f in [flux, flux_err]:
            assert np.sum(~np.isnan(flux)) == xtrace.size
        # test npix is zero there too
        assert np.count_nonzero(npix) == xtrace.size

        # test that most of npix are equal to width (Not all, because of NaN mask, but NaN fraction)
        # is small enough that it should still be the most represented count for such a small width
        unique, counts = np.unique(npix, return_counts=True)
        assert np.isclose(unique[np.argmax(counts)], WIDTH)


def test_estim_error_nearest_data(imagemodel, mask_trace_profile):
    data, err = imagemodel

    for order in [0, 1]:
        # this has bad pixels set to 1, ONLY within the spectral trace.
        # everything else is zero, i.e., regions outside trace and good data
        pix_to_estim = np.zeros(data.shape, dtype="bool")
        pix_to_estim[np.isnan(data)] = 1
        pix_to_estim[mask_trace_profile[order] == 1] = 0

        # this has bad pixels set to 0, and regions outside trace set to 0, and good data 1
        valid_pix = ~mask_trace_profile[order]
        valid_pix[pix_to_estim] = 0

        err_out = estim_error_nearest_data(err, data, pix_to_estim, valid_pix)

        # test that all replaced values are positive and no NaNs are left
        assert np.sum(np.isnan(err_out)) == 0
        assert np.all(err_out > 0)

        # # test that the replaced pixels are not statistical outliers c.f. the other pixels
        # replaced_pix = err_out[pix_to_estim]
        # original_pix = err_out[valid_pix]
        # diff = np.mean(replaced_pix)/np.mean(original_pix)
        # assert np.isclose(diff, 1, rtol=0.5) # assert False
        # # this fails because estim_error_nearest_data takes the smaller of the two nearest
        # # good pixels, which leads to the replaced pixels being lower on average
        # # than the original pixels.
