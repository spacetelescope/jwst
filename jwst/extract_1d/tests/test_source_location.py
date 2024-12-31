import logging
import numpy as np
import pytest

from jwst.extract_1d import source_location as ex
from jwst.tests.helpers import LogWatcher


@pytest.fixture
def log_watcher(monkeypatch):
    # Set a log watcher to check for a log message at any level
    # in the extract_1d.extract module
    watcher = LogWatcher('')
    logger = logging.getLogger('jwst.extract_1d.source_location')
    for level in ['debug', 'info', 'warning', 'error']:
        monkeypatch.setattr(logger, level, watcher)
    return watcher


@pytest.mark.parametrize('resampled', [True, False])
@pytest.mark.parametrize('is_slit', [True, False])
@pytest.mark.parametrize('missing_bbox', [True, False])
def test_location_from_wcs_nirspec(
        monkeypatch, mock_nirspec_fs_one_slit, resampled, is_slit, missing_bbox):
    model = mock_nirspec_fs_one_slit

    if not resampled:
        # mock available frames, so it looks like unresampled cal data
        monkeypatch.setattr(model.meta.wcs, 'available_frames', ['gwa'])

    if missing_bbox:
        # mock a missing bounding box - should have same results
        # for the test data
        monkeypatch.setattr(model.meta.wcs, 'bounding_box', None)

    if is_slit:
        middle, middle_wl, location, trace = ex.location_from_wcs(model, model)
    else:
        middle, middle_wl, location, trace = ex.location_from_wcs(model, None)

    # middle pixel is center of dispersion axis
    assert middle == int((model.data.shape[1] - 1) / 2)

    # middle wavelength is the wavelength at that point, from the mock wcs
    assert np.isclose(middle_wl, 7.745)

    # location is 1.0 - from the mocked transform function
    assert location == 1.0

    # trace is the same, in an array
    assert np.all(trace == 1.0)


@pytest.mark.parametrize('is_slit', [True, False])
def test_location_from_wcs_miri(monkeypatch, mock_miri_lrs_fs, is_slit):
    model = mock_miri_lrs_fs

    # monkey patch in a transform for the wcs
    def radec2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 1.0, 0.0
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'backward_transform', radec2det())

    # mock the trace function
    def mock_trace(*args, **kwargs):
        return np.full(model.data.shape[-2], 1.0)

    monkeypatch.setattr(ex, '_miri_trace_from_wcs', mock_trace)

    # Get the slit center from the WCS
    if is_slit:
        middle, middle_wl, location, trace = ex.location_from_wcs(model, model)
    else:
        middle, middle_wl, location, trace = ex.location_from_wcs(model, None)

    # middle pixel is center of dispersion axis
    assert middle == int((model.data.shape[0] - 1) / 2)

    # middle wavelength is the wavelength at that point, from the mock wcs
    assert np.isclose(middle_wl, 7.255)

    # location is 1.0 - from the mocked transform function
    assert location == 1.0

    # trace is the same, in an array
    assert np.all(trace == 1.0)


def test_location_from_wcs_missing_data(mock_miri_lrs_fs, log_watcher):
    model = mock_miri_lrs_fs
    model.meta.wcs.backward_transform = None

    # model is missing WCS information - None values are returned
    log_watcher.message = "Dithered pointing location not found"
    result = ex.location_from_wcs(model, None)
    assert result == (None, None, None, None)
    log_watcher.assert_seen()


def test_location_from_wcs_wrong_exptype(mock_niriss_soss, log_watcher):
    # model is not a handled exposure type
    log_watcher.message = "Source position cannot be found for EXP_TYPE"
    result = ex.location_from_wcs(mock_niriss_soss, None)
    assert result == (None, None, None, None)
    log_watcher.assert_seen()


def test_location_from_wcs_bad_location(
        monkeypatch, mock_nirspec_fs_one_slit, log_watcher):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, np.nan
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'get_transform', slit2det)

    # WCS transform returns NaN for the location
    log_watcher.message = "Source position could not be determined"
    result = ex.location_from_wcs(model, None)
    assert result == (None, None, None, None)
    log_watcher.assert_seen()


def test_location_from_wcs_location_out_of_range(
        monkeypatch, mock_nirspec_fs_one_slit, log_watcher):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, 2000
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'get_transform', slit2det)

    # mock the trace function
    def mock_trace(*args, **kwargs):
        return np.full(model.data.shape[-1], 1.0)

    monkeypatch.setattr(ex, '_nirspec_trace_from_wcs', mock_trace)

    # WCS transform a value outside the bounding box
    log_watcher.message = "outside the bounding box"
    result = ex.location_from_wcs(model, None)
    assert result == (None, None, None, None)
    log_watcher.assert_seen()


def test_nirspec_trace_from_wcs(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    trace = ex._nirspec_trace_from_wcs(model.data.shape, model.meta.wcs.bounding_box,
                                       model.meta.wcs, 1.0, 1.0)
    # mocked model contains some mock transforms as well - all ones are expected
    assert np.all(trace == np.ones(model.data.shape[-1]))


def test_miri_trace_from_wcs(mock_miri_lrs_fs):
    model = mock_miri_lrs_fs
    trace = ex._miri_trace_from_wcs(model.data.shape, model.meta.wcs.bounding_box,
                                    model.meta.wcs, 1.0, 1.0)

    # mocked model contains some mock transforms as well - all ones are expected
    assert np.all(trace == np.ones(model.data.shape[-1]))
