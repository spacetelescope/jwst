import numpy as np
import pytest

from jwst.extract_1d import source_location as sl


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_middle_from_wcs_constant_wl(dispaxis):
    # mock a wcs that returns a constant wavelength
    def mock_wcs(x, y):
        return None, None, np.full(x.shape, 10.0)

    bbox = (-0.5, 9.5), (-0.5, 9.5)

    md, mx, mw = sl.middle_from_wcs(mock_wcs, bbox, dispaxis)

    # middle for dispersion and cross-dispersion are at the center, 4.5
    assert md == 4.5
    assert mx == 4.5

    # middle wavelength is the constant value
    assert mw == 10.0


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_middle_from_wcs_variable_wl(dispaxis):
    # mock a wcs that returns a variable wavelength
    def mock_wcs(x, y):
        return None, None, np.arange(x.size, dtype=float).reshape(x.shape)

    bbox = (-0.5, 9.5), (-0.5, 9.5)

    md, mx, mw = sl.middle_from_wcs(mock_wcs, bbox, dispaxis)

    # middle for dispersion, cross-dispersion, and wavelength are all 4.5
    assert md == 4.5
    assert mx == 4.5
    assert mw == 4.5


@pytest.mark.parametrize("resampled", [True, False])
@pytest.mark.parametrize("is_slit", [True, False])
@pytest.mark.parametrize("missing_bbox", [True, False])
def test_location_from_wcs_nirspec(
    monkeypatch, mock_nirspec_fs_one_slit, resampled, is_slit, missing_bbox
):
    model = mock_nirspec_fs_one_slit

    if not resampled:
        # mock available frames, so it looks like unresampled cal data
        monkeypatch.setattr(model.meta.wcs, "available_frames", ["gwa"])

    if missing_bbox:
        # mock a missing bounding box - should have same results
        # for the test data
        monkeypatch.setattr(model.meta.wcs, "bounding_box", None)

    if is_slit:
        middle, middle_wl, location, trace = sl.location_from_wcs(model, model)
    else:
        middle, middle_wl, location, trace = sl.location_from_wcs(model, None)

    # middle pixel is center of dispersion axis
    assert middle == int((model.data.shape[1] - 1) / 2)

    # middle wavelength is the wavelength at that point, from the mock wcs
    assert np.isclose(middle_wl, 7.745)

    # location is 1.0 - from the mocked transform function
    assert location == 1.0

    # trace is the same, in an array
    assert np.all(trace == 1.0)


@pytest.mark.parametrize("is_slit", [True, False])
def test_location_from_wcs_miri(monkeypatch, mock_miri_lrs_fs, is_slit):
    model = mock_miri_lrs_fs

    # monkey patch in a transform for the wcs
    def radec2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 1.0, 0.0

        return return_one

    monkeypatch.setattr(model.meta.wcs, "backward_transform", radec2det())

    # mock the trace function
    def mock_trace(*args, **kwargs):
        return np.full(model.data.shape[-2], 1.0)

    monkeypatch.setattr(sl, "_miri_trace_from_wcs", mock_trace)

    # Get the slit center from the WCS
    if is_slit:
        middle, middle_wl, location, trace = sl.location_from_wcs(model, model)
    else:
        middle, middle_wl, location, trace = sl.location_from_wcs(model, None)

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
    watcher = log_watcher(
        "jwst.extract_1d.source_location",
        message="Dithered pointing location not found",
        level="warning",
    )

    # model is missing WCS information - None values are returned
    result = sl.location_from_wcs(model, None)
    assert result == (None, None, None, None)
    watcher.assert_seen()


def test_location_from_wcs_wrong_exptype(mock_niriss_soss, log_watcher):
    # model is not a handled exposure type
    watcher = log_watcher(
        "jwst.extract_1d.source_location",
        message="Source position cannot be found for EXP_TYPE",
        level="warning",
    )
    result = sl.location_from_wcs(mock_niriss_soss, None)
    assert result == (None, None, None, None)
    watcher.assert_seen()


def test_location_from_wcs_bad_location(monkeypatch, mock_nirspec_fs_one_slit, log_watcher):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, np.nan

        return return_one

    monkeypatch.setattr(model.meta.wcs, "get_transform", slit2det)

    # WCS transform returns NaN for the location
    watcher = log_watcher(
        "jwst.extract_1d.source_location",
        message="Source position could not be determined",
        level="warning",
    )
    result = sl.location_from_wcs(model, None)
    assert result == (None, None, None, None)
    watcher.assert_seen()


def test_location_from_wcs_location_out_of_range(
    monkeypatch, mock_nirspec_fs_one_slit, log_watcher
):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, 2000

        return return_one

    monkeypatch.setattr(model.meta.wcs, "get_transform", slit2det)

    # mock the trace function
    def mock_trace(*args, **kwargs):
        return np.full(model.data.shape[-1], 1.0)

    monkeypatch.setattr(sl, "_nirspec_trace_from_wcs", mock_trace)

    # WCS transform a value outside the bounding box
    watcher = log_watcher(
        "jwst.extract_1d.source_location", message="outside the bounding box", level="warning"
    )
    result = sl.location_from_wcs(model, None)
    assert result == (None, None, None, None)
    watcher.assert_seen()


def test_nirspec_trace_from_wcs(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    trace = sl._nirspec_trace_from_wcs(
        model.data.shape, model.meta.wcs.bounding_box, model.meta.wcs, 1.0, 1.0
    )
    # mocked model contains some mock transforms as well - all ones are expected
    assert np.all(trace == np.ones(model.data.shape[1]))


def test_miri_trace_from_wcs(mock_miri_lrs_fs):
    model = mock_miri_lrs_fs
    trace = sl._miri_trace_from_wcs(
        model.data.shape, model.meta.wcs.bounding_box, model.meta.wcs, 1.0, 1.0
    )

    # mocked model contains some mock transforms as well - all ones are expected
    assert np.all(trace == np.ones(model.data.shape[0]))


def test_trace_from_wcs_nirspec(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    trace = sl.trace_from_wcs(
        "NRS_FIXEDSLIT", model.data.shape, model.meta.wcs.bounding_box, model.meta.wcs, 1.0, 1.0, 1
    )

    # mocked model contains some mock transforms as well - all ones are expected
    assert np.all(trace == np.ones(model.data.shape[1]))


def test_trace_from_wcs_miri(mock_miri_lrs_fs):
    model = mock_miri_lrs_fs
    trace = sl.trace_from_wcs(
        "MIR_LRS-FIXEDSLIT",
        model.data.shape,
        model.meta.wcs.bounding_box,
        model.meta.wcs,
        1.0,
        1.0,
        2,
    )

    # mocked model contains some mock transforms as well - all ones are expected
    assert np.all(trace == np.ones(model.data.shape[0]))


def test_trace_from_wcs_other_horizontal():
    exp_type = "ANY"
    shape = (10, 20)
    bbox = (1.5, 17.5), (1.5, 8.5)
    wcs = None
    source_x = 2.0
    source_y = 4.0
    dispaxis = 1

    trace = sl.trace_from_wcs(exp_type, shape, bbox, wcs, source_x, source_y, dispaxis)

    # trace matches dispersion dimension (x)
    assert trace.shape == (shape[1],)

    # trace is full of the source_y position, except outside the bounding box,
    # where it is NaN
    assert np.all(trace[2:18] == source_y)
    assert np.all(np.isnan(trace[:2]))
    assert np.all(np.isnan(trace[18:]))


def test_trace_from_wcs_other_vertical():
    exp_type = "ANY"
    shape = (10, 20)
    bbox = (1.5, 17.5), (1.5, 8.5)
    wcs = None
    source_x = 2.0
    source_y = 4.0
    dispaxis = 2

    trace = sl.trace_from_wcs(exp_type, shape, bbox, wcs, source_x, source_y, dispaxis)

    # trace matches dispersion dimension (y)
    assert trace.shape == (shape[0],)

    # trace is full of the source_x position, except outside the bounding box,
    # where it is NaN
    assert np.all(trace[2:9] == source_x)
    assert np.all(np.isnan(trace[:2]))
    assert np.all(np.isnan(trace[9:]))


def test_nod_pair_location_nirspec(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    middle_wl = 7.5

    nod_center = sl.nod_pair_location(model, middle_wl)

    # for mock transforms, 1.0 is expected
    assert nod_center == 1.0


def test_nod_pair_location_nirspec_unresampled(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    middle_wl = 7.5
    model.meta.wcs.available_frames = ["gwa"]

    nod_center = sl.nod_pair_location(model, middle_wl)

    # for mock transforms, 1.0 is expected
    assert nod_center == 1.0


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_nod_pair_location_miri(mock_miri_lrs_fs, dispaxis):
    model = mock_miri_lrs_fs
    middle_wl = 7.5

    nod_center = sl.nod_pair_location(model, middle_wl)

    # for mock transforms as is, NaN is expected
    assert np.isnan(nod_center)

    # mock v2v3 transform
    model.meta.wcs.available_frames = ["v2v3"]
    model.meta.wcsinfo.v3yangle = 1.0
    model.meta.wcsinfo.v2_ref = 1.0
    model.meta.wcsinfo.v3_ref = 1.0
    model.meta.wcsinfo.vparity = 1
    model.meta.dither.x_offset = 1.0
    model.meta.dither.y_offset = 1.0
    model.meta.wcsinfo.dispersion_direction = dispaxis
    nod_center = sl.nod_pair_location(model, middle_wl)

    # the final backward transform mock returns (1.0, 0.0),
    # so the location reported will be x=1.0 if vertical, y=0.0 if horizontal
    if dispaxis == 2:
        assert nod_center == 1.0
    else:
        assert nod_center == 0.0
