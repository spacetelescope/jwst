import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.stats import sigma_clipped_stats

from jwst.wfss_contam.observations import background_subtract, _select_ids, Observation
import stdatamodels.jwst.datamodels as dm


@pytest.fixture(scope="module")
def observation(direct_image_with_gradient, segmentation_map, grism_wcs):
    """
    set up observation object with mock data.
    direct_image_with_gradient still needs to be run to produce the file,
    even though it is not called directly
    """
    filter_name = "F200W"
    seg = segmentation_map.data
    all_ids = np.array(list(set(np.ravel(seg))))
    source_ids = all_ids[50:52]
    obs = Observation(
        direct_image_with_gradient.data,
        segmentation_map,
        grism_wcs,
        filter_name,
        source_id=source_ids,
    )
    return obs


@pytest.mark.parametrize(
    "source_id, expected", [(None, [1, 2, 3]), (2, [2]), (np.array([1, 3]), [1, 3])]
)
def test_select_ids(source_id, expected):
    all_ids = [1, 2, 3]
    assert _select_ids(source_id, all_ids) == expected
    assert isinstance(_select_ids(source_id, all_ids), list)


def test_select_ids_expected_raises():
    with pytest.raises(ValueError):
        _select_ids("all", [1, 2, 3])


def test_background_subtract(direct_image_with_gradient):
    data = direct_image_with_gradient.data
    subtracted_data = background_subtract(data)
    mean, _median, stddev = sigma_clipped_stats(subtracted_data, sigma=3.0)
    assert_allclose(mean, 0.0, atol=0.2 * stddev)


def test_disperse_order(observation):
    obs = observation
    order = 1
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)

    # manually change x,y offset because took transform from a real direct image, with different
    # pixel 0,0 than the mock data. This puts i=1, order 1 onto the real grism image
    obs.xoffset = 2200
    obs.yoffset = 1000

    # shorten pixel list to make this test take less time
    obs.disperse_order(order, wmin, wmax, sens_waves, sens_resp)

    # test simulated image. should be majority but not all zeros
    assert obs.simulated_image.shape == obs.dims
    assert not np.allclose(obs.simulated_image, 0.0)
    assert np.median(obs.simulated_image) == 0.0

    # test simulated slits and their associated metadata
    # only the second of the two obs ids is in the simulated image
    assert type(obs.simulated_slits) == dm.MultiSlitModel
    assert len(obs.simulated_slits.slits) == 1
    slit = obs.simulated_slits.slits[0]
    # check metadata
    assert slit.name == "source_51"
    assert slit.data.shape == (slit.ysize, slit.xsize)

    # check for regression by hard-coding one value of slit.data
    assert np.isclose(slit.data[10,10], 111.29376)
