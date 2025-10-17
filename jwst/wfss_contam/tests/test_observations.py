import copy
import logging

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.stats import sigma_clipped_stats
from numpy.testing import assert_allclose

from jwst.tests.helpers import LogWatcher
from jwst.wfss_contam.observations import Observation, _select_ids, background_subtract


@pytest.fixture
def observation(direct_image_with_gradient, segmentation_map, grism_wcs):
    """
    set up observation object with mock data.
    direct_image_with_gradient still needs to be run to produce the file,
    even though it is not called directly
    """
    seg = segmentation_map.data
    return Observation(
        direct_image_with_gradient.data,
        segmentation_map.data,
        grism_wcs,
        direct_image_with_gradient.meta.wcs,
        phot_per_lam=False,
    )


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


def test_chunk_sources(observation, segmentation_map, monkeypatch):
    obs = copy.deepcopy(observation)
    order = 1
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_response = np.ones(100)

    seg = segmentation_map.data
    all_ids = np.array(list(set(np.ravel(seg))))
    source_ids = all_ids[50:60]

    # find largest source out of source_ids, then make max_pixels smaller than that
    source_ids_per_pixel = obs.source_ids_per_pixel[np.isin(obs.source_ids_per_pixel, source_ids)]
    ids, n_pix_per_sources = np.unique(source_ids_per_pixel, return_counts=True)
    max_pixels = np.max(n_pix_per_sources) - 1  # to trigger the warning
    bad_id = ids[n_pix_per_sources > max_pixels][0]

    # ensure warning is emitted for source that is too large
    watcher = LogWatcher(
        f"Source {bad_id} has {np.max(n_pix_per_sources)} pixels, which exceeds the maximum"
    )
    monkeypatch.setattr(logging.getLogger("jwst.wfss_contam.observations"), "warning", watcher)
    disperse_args = obs.chunk_sources(
        order,
        wmin,
        wmax,
        sens_waves,
        sens_response,
        selected_ids=source_ids,
        max_pixels=max_pixels,
    )
    watcher.assert_seen()

    # two of the sources were too large and skipped
    assert len(disperse_args) == 8


def test_disperse_order(observation, segmentation_map):
    obs = copy.deepcopy(observation)
    order = 1
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones_like(sens_waves)

    seg = segmentation_map.data
    all_ids = np.array(list(set(np.ravel(seg))))
    source_ids = all_ids[50:60]

    # shorten pixel list to make this test take less time
    obs.disperse_order(order, wmin, wmax, sens_waves, sens_resp, selected_ids=source_ids)

    # test simulated image. should be majority but not all zeros
    assert obs.simulated_image.shape == obs.dims
    assert not np.allclose(obs.simulated_image, 0.0)
    assert np.median(obs.simulated_image) == 0.0

    # test simulated slits and their associated metadata
    assert type(obs.simulated_slits) == dm.MultiSlitModel
    # two of the slits fall off the detector and do not end up here
    assert len(obs.simulated_slits.slits) == 8
    slit = obs.simulated_slits.slits[1]
    # check metadata
    assert slit.name == "source_51"
    assert slit.data.shape == (slit.ysize, slit.xsize)

    # check for regression by hard-coding one value of slit.data
    assert np.isclose(slit.data[5, 60], 20.996877670288086)
