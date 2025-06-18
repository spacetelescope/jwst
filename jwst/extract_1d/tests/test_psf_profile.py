import numpy as np
import pytest
from stdatamodels.jwst.datamodels import SpecPsfModel

from jwst.extract_1d import psf_profile as pp


@pytest.mark.parametrize("exp_type", ["MIR_LRS-FIXEDSLIT", "NRS_FIXEDSLIT", "UNKNOWN"])
def test_open_psf(psf_reference_file, exp_type):
    # for any exptype, a model that can be read
    # as SpecPsfModel will be, since it's the only
    # one implemented so far
    with pp.open_psf(psf_reference_file, exp_type=exp_type) as model:
        assert isinstance(model, SpecPsfModel)


def test_open_psf_fail():
    with pytest.raises(NotImplementedError, match="could not be read"):
        pp.open_psf("bad_file", "UNKNOWN")


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_normalize_profile(nod_profile, dispaxis):
    profile = 2 * nod_profile
    if dispaxis == 2:
        profile = profile.T
    pp._normalize_profile(profile, dispaxis)
    assert np.allclose(np.sum(profile, axis=dispaxis - 1), 1.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_normalize_profile_with_nans(nod_profile, dispaxis):
    profile = -1 * nod_profile
    profile[10, :] = np.nan
    if dispaxis == 2:
        profile = profile.T

    pp._normalize_profile(profile, dispaxis)
    assert np.allclose(np.sum(profile, axis=dispaxis - 1), 1.0)
    assert np.all(np.isfinite(profile))


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_make_cutout_profile_default(psf_reference, dispaxis):
    data_shape = psf_reference.data.shape
    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]

    psf_subpix = psf_reference.meta.psf.subpix
    profiles = pp._make_cutout_profile(xidx, yidx, psf_subpix, psf_reference.data, dispaxis)
    assert len(profiles) == 1
    assert profiles[0].shape == data_shape

    # No shift, profile is uniform and normalized to cross-dispersion size
    if dispaxis == 1:
        assert np.all(profiles[0] == 1 / data_shape[0])
    else:
        assert np.all(profiles[0] == 1 / data_shape[1])


@pytest.mark.parametrize("dispaxis", [1, 2])
@pytest.mark.parametrize("extra_shift", [1, 2])
def test_make_cutout_profile_shift_down(psf_reference, dispaxis, extra_shift):
    data_shape = psf_reference.data.shape
    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    psf_subpix = psf_reference.meta.psf.subpix

    profiles = pp._make_cutout_profile(
        xidx, yidx, psf_subpix, psf_reference.data, dispaxis, extra_shift=extra_shift
    )
    assert len(profiles) == 1
    assert profiles[0].shape == data_shape

    # Profile is shifted down by shift amount
    if dispaxis == 1:
        nn = data_shape[0] - extra_shift
        assert np.all(profiles[0][:-extra_shift, :] == 1 / nn)
        assert np.all(profiles[0][-extra_shift:, :] == 0.0)
    else:
        nn = data_shape[1] - extra_shift
        assert np.all(profiles[0][:, :-extra_shift] == 1 / nn)
        assert np.all(profiles[0][:, -extra_shift:] == 0.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
@pytest.mark.parametrize("extra_shift", [-1, -2])
def test_make_cutout_profile_shift_up(psf_reference, dispaxis, extra_shift):
    data_shape = psf_reference.data.shape
    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    psf_subpix = psf_reference.meta.psf.subpix

    profiles = pp._make_cutout_profile(
        xidx, yidx, psf_subpix, psf_reference.data, dispaxis, extra_shift=extra_shift
    )
    assert len(profiles) == 1
    assert profiles[0].shape == data_shape

    # Profile is shifted up by shift amount
    if dispaxis == 1:
        nn = data_shape[0] + extra_shift
        assert np.all(profiles[0][-extra_shift:, :] == 1 / nn)
        assert np.all(profiles[0][:-extra_shift, :] == 0.0)
    else:
        nn = data_shape[1] + extra_shift
        assert np.all(profiles[0][:, -extra_shift:] == 1 / nn)
        assert np.all(profiles[0][:, :-extra_shift] == 0.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_make_cutout_profile_with_nod(psf_reference, dispaxis):
    data_shape = psf_reference.data.shape
    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    psf_subpix = psf_reference.meta.psf.subpix

    offset = 2
    profiles = pp._make_cutout_profile(
        xidx, yidx, psf_subpix, psf_reference.data, dispaxis, nod_offset=offset
    )
    assert len(profiles) == 2
    source, nod = profiles
    assert source.shape == data_shape
    assert nod.shape == data_shape

    # Profile is uniform, nod profile is the same, but shifted
    # down by offset and multiplied by -1
    if dispaxis == 1:
        assert np.all(source == 1 / data_shape[0])

        nn = data_shape[0] - offset
        assert np.all(nod[:-offset, :] == -1 / nn)
        assert np.all(nod[-offset:, :] == 0.0)
    else:
        assert np.all(source == 1 / data_shape[1])

        nn = data_shape[1] - offset
        assert np.all(nod[:, :-offset] == -1 / nn)
        assert np.all(nod[:, -offset:] == 0.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_profile_residual(psf_reference, dispaxis):
    data_shape = (50, 50)
    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    psf_subpix = psf_reference.meta.psf.subpix

    # Set data to all ones, so residual should be zero
    # when background is not fit
    data = np.full(data_shape, 1.0)
    var = np.full(data_shape, 0.01)

    param = [0, None]
    residual = pp._profile_residual(
        param, data, var, xidx, yidx, psf_subpix, psf_reference.data, dispaxis, fit_bkg=False
    )
    assert np.isclose(residual, 0.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_profile_residual_with_bkg(psf_reference, dispaxis):
    data_shape = (50, 50)
    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    psf_subpix = psf_reference.meta.psf.subpix

    # Set data to all ones, so it is all background - residual
    # should be all of the data
    data = np.full(data_shape, 1.0)
    var = np.full(data_shape, 0.01)

    param = [0, None]
    residual = pp._profile_residual(
        param, data, var, xidx, yidx, psf_subpix, psf_reference.data, dispaxis, fit_bkg=True
    )
    assert np.isclose(residual, np.sum(data**2 / var))


@pytest.mark.parametrize("use_trace", [True, False])
def test_psf_profile(mock_miri_lrs_fs, psf_reference_file, use_trace):
    data_shape = mock_miri_lrs_fs.data.shape
    if use_trace:
        # Centered trace - should match None behavior
        trace = np.full(data_shape[0], (data_shape[1] - 1) / 2.0)
    else:
        # Avoid miri specific behavior for trace - the mock WCS is not complete enough
        mock_miri_lrs_fs.meta.exposure.type = "ANY"
        trace = None

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = mock_miri_lrs_fs.meta.wcs(xidx, yidx)

    profiles, lower, upper = pp.psf_profile(
        mock_miri_lrs_fs,
        trace,
        wl_array,
        psf_reference_file,
        optimize_shifts=False,
        model_nod_pair=False,
    )

    # no nod profile, data cutout matches full shape
    assert len(profiles) == 1
    assert profiles[0].shape == data_shape
    assert lower == 0
    assert upper == data_shape[1]

    # profile is uniform and centered
    assert np.allclose(profiles[0], 1 / data_shape[1])


@pytest.mark.parametrize("use_trace", [True, False])
def test_psf_profile_multi_int(mock_nirspec_bots, psf_reference_file, use_trace):
    data_shape = mock_nirspec_bots.data.shape[-2:]

    if use_trace:
        trace = np.full(data_shape[1], (data_shape[0] - 1) / 2.0)
    else:
        # Avoid nirspec specific behavior for trace -
        # the mock WCS is not complete enough
        mock_nirspec_bots.meta.exposure.type = "ANY"
        trace = None

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = mock_nirspec_bots.meta.wcs(xidx, yidx)

    profiles, lower, upper = pp.psf_profile(
        mock_nirspec_bots,
        trace,
        wl_array,
        psf_reference_file,
        optimize_shifts=False,
        model_nod_pair=False,
    )

    # no nod profile, data cutout matches full shape
    assert len(profiles) == 1
    assert profiles[0].shape == data_shape
    assert lower == 0
    assert upper == data_shape[0]

    # profile is uniform and centered
    assert np.allclose(profiles[0], 1 / data_shape[1])


def test_psf_profile_model_nod(monkeypatch, mock_miri_lrs_fs, psf_reference_file):
    model = mock_miri_lrs_fs
    data_shape = model.data.shape
    trace = np.full(data_shape[0], (data_shape[1] - 1) / 2.0)

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    # mock nod subtraction
    model.meta.cal_step.bkg_subtract = "COMPLETE"
    model.meta.dither.primary_type = "2-POINT-NOD"

    # mock a nod position at the opposite end of the array
    def mock_nod(*args, **kwargs):
        return 48.0

    monkeypatch.setattr(pp, "nod_pair_location", mock_nod)

    profiles, lower, upper = pp.psf_profile(
        model, trace, wl_array, psf_reference_file, optimize_shifts=False, model_nod_pair=True
    )

    # now returns nod profile
    assert len(profiles) == 2
    source, nod = profiles
    assert source.shape == data_shape
    assert nod.shape == data_shape

    # source profile is uniform and centered
    assert np.allclose(source, 1 / data_shape[1])

    # nod profile is centered at the end of the array and has a negative value
    assert np.allclose(nod[:, -2:], -1 / 26)


def test_psf_profile_model_nod_no_trace(mock_miri_lrs_fs, psf_reference_file, log_watcher):
    model = mock_miri_lrs_fs
    model.meta.exposure.type = "ANY"
    data_shape = model.data.shape
    trace = None

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    watcher = log_watcher(
        "jwst.extract_1d.psf_profile", message="Cannot model a negative nod without position"
    )
    profiles, lower, upper = pp.psf_profile(
        model, trace, wl_array, psf_reference_file, optimize_shifts=False, model_nod_pair=True
    )
    watcher.assert_seen()

    # does not return nod profile
    assert len(profiles) == 1


def test_psf_profile_model_nod_not_subtracted(mock_miri_lrs_fs, psf_reference_file, log_watcher):
    model = mock_miri_lrs_fs
    data_shape = model.data.shape
    trace = np.full(data_shape[0], (data_shape[1] - 1) / 2.0)

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    watcher = log_watcher("jwst.extract_1d.psf_profile", message="data was not nod-subtracted")
    profiles, lower, upper = pp.psf_profile(
        model, trace, wl_array, psf_reference_file, optimize_shifts=False, model_nod_pair=True
    )
    watcher.assert_seen()

    # does not return nod profile
    assert len(profiles) == 1


def test_psf_profile_model_nod_wrong_pattern(mock_miri_lrs_fs, psf_reference_file, log_watcher):
    model = mock_miri_lrs_fs
    data_shape = model.data.shape
    trace = np.full(data_shape[0], (data_shape[1] - 1) / 2.0)
    model.meta.cal_step.bkg_subtract = "COMPLETE"

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    watcher = log_watcher("jwst.extract_1d.psf_profile", message="data was not a two-point nod")
    profiles, lower, upper = pp.psf_profile(
        model, trace, wl_array, psf_reference_file, optimize_shifts=False, model_nod_pair=True
    )
    watcher.assert_seen()

    # does not return nod profile
    assert len(profiles) == 1


def test_psf_profile_model_nod_bad_position(mock_miri_lrs_fs, psf_reference_file, log_watcher):
    model = mock_miri_lrs_fs
    data_shape = model.data.shape
    trace = np.full(data_shape[0], (data_shape[1] - 1) / 2.0)
    model.meta.cal_step.bkg_subtract = "COMPLETE"
    model.meta.dither.primary_type = "2-POINT-NOD"

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    watcher = log_watcher(
        "jwst.extract_1d.psf_profile", message="Nod center could not be estimated"
    )
    profiles, lower, upper = pp.psf_profile(
        model, trace, wl_array, psf_reference_file, optimize_shifts=False, model_nod_pair=True
    )
    watcher.assert_seen()

    # does not return nod profile
    assert len(profiles) == 1


@pytest.mark.parametrize("start_offset", [0.0, 0.1, 2.0, -1.0])
def test_psf_profile_optimize(
    mock_miri_lrs_fs, psf_reference_file_with_source, start_offset, log_watcher
):
    model = mock_miri_lrs_fs
    data_shape = model.data.shape
    trace = np.full(data_shape[0], (data_shape[1] - 1) / 2.0)
    trace += start_offset

    # add a peak to the center of the data, matching the model
    model.data[:] = 0.01
    model.data[:, 24:27] += 1.0
    model.var_rnoise[:] = 0.1

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    watcher = log_watcher(
        "jwst.extract_1d.psf_profile", message="Centering profile on spectrum at 24.5"
    )
    profiles, lower, upper = pp.psf_profile(
        model,
        trace,
        wl_array,
        psf_reference_file_with_source,
        optimize_shifts=True,
        model_nod_pair=False,
    )
    watcher.assert_seen()

    # profile is centered at 24.5
    profile = profiles[0]
    assert np.allclose(profile[:, 24:27], 1 / 3, atol=1e-4)
    assert np.allclose(profile[:, :24], 0.0, atol=1e-4)
    assert np.allclose(profile[:, 27:], 0.0, atol=1e-4)


@pytest.mark.parametrize("start_offset", [0.0, 0.1, 2.0, -1.0])
def test_psf_profile_optimize_with_nod(
    monkeypatch, mock_miri_lrs_fs, psf_reference_file_with_source, start_offset, log_watcher
):
    model = mock_miri_lrs_fs
    data_shape = model.data.shape

    # mock nod subtraction
    model.meta.cal_step.bkg_subtract = "COMPLETE"
    model.meta.dither.primary_type = "2-POINT-NOD"

    # trace at pixel 9.5
    trace = np.full(data_shape[0], 9.5)

    # with an offset to optimize out
    trace += start_offset

    # mock a nod position at the opposite end of the array,
    # with an extra offset from truth, and also from the trace
    def mock_nod(*args, **kwargs):
        return 39.5 + start_offset + 0.1

    monkeypatch.setattr(pp, "nod_pair_location", mock_nod)

    # add a peak at pixel 10, negative peak at pixel 40
    model.data[:] = 0.0
    model.data[:, 9:12] += 1.0
    model.data[:, 39:42] -= 1.0
    model.var_rnoise[:] = 0.1

    yidx, xidx = np.mgrid[: data_shape[0], : data_shape[1]]
    _, _, wl_array = model.meta.wcs(xidx, yidx)

    watcher = log_watcher(
        "jwst.extract_1d.psf_profile", message="Also modeling a negative trace at 39.50"
    )
    profiles, lower, upper = pp.psf_profile(
        model,
        trace,
        wl_array,
        psf_reference_file_with_source,
        optimize_shifts=True,
        model_nod_pair=True,
    )
    watcher.assert_seen()

    # profile is centered at 10
    source, nod = profiles
    assert np.allclose(source[:, 9:12], 1 / 3, atol=1e-4)
    assert np.allclose(source[:, :9], 0.0, atol=1e-4)
    assert np.allclose(source[:, 12:], 0.0, atol=1e-4)

    # nod is centered at 40
    assert np.allclose(nod[:, 39:42], -1 / 3, atol=1e-4)
    assert np.allclose(nod[:, :39], 0.0, atol=1e-4)
    assert np.allclose(nod[:, 42:], 0.0, atol=1e-4)
