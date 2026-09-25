import os
from glob import glob

import numpy as np
import pytest
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.dqflags import pixel as flags

from jwst.datamodels import ModelContainer
from jwst.pixel_replace.pixel_replace import PixelReplacement
from jwst.pixel_replace.pixel_replace_step import PixelReplaceStep
from jwst.pixel_replace.tests import helpers


@pytest.fixture(scope="module")
def nirspec_tso():
    model, bad_idx = helpers.nirspec_tso()
    yield model, bad_idx
    model.close()


@pytest.fixture(scope="module")
def nirspec_fs_slitmodel():
    model, bad_idx = helpers.nirspec_fs_slitmodel()
    yield model, bad_idx
    model.close()


@pytest.fixture(scope="module")
def miri_lrs():
    model, bad_idx = helpers.miri_lrs()
    yield model, bad_idx
    model.close()


@pytest.fixture(scope="module")
def miri_mrs():
    model, bad_idx = helpers.miri_mrs()
    yield model, bad_idx
    model.close()


@pytest.fixture(scope="module")
def nirspec_ifu():
    model, bad_idx = helpers.nirspec_ifu()
    yield model, bad_idx
    model.close()


@pytest.fixture(scope="module")
def nirspec_msa_multislit():
    model, bad_idx = helpers.nirspec_msa_multislit()
    yield model, bad_idx
    model.close()


@pytest.mark.parametrize(
    "dataset", ["nirspec_tso", "nirspec_fs_slitmodel", "miri_lrs", "miri_mrs", "nirspec_ifu"]
)
@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_pixel_replace_no_container(request, dataset, algorithm):
    """
    Test pixel replace for modes with no container.

    This includes ImageModel, SlitModel, and IFUImageModel.
    """
    input_model, bad_idx = request.getfixturevalue(dataset)

    # for this simple case, the results from either algorithm should
    # be the same
    result = PixelReplaceStep.call(input_model, skip=False, algorithm=algorithm)

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result, ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result, ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model, ext)[bad_idx])
        assert getattr(result, ext)[bad_idx] == 1.0

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.dq[bad_idx] == (
        input_model.dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.dq[..., :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.dq[..., 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    # Step is recorded as complete
    assert result.meta.cal_step.pixel_replace == "COMPLETE"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.pixel_replace is None

    result.close()
    input_model.close()


@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_pixel_replace_multislit(nirspec_msa_multislit, algorithm):
    """Test pixel replace for multislit modes."""
    input_model, bad_idx = nirspec_msa_multislit

    # for this simple case, the results from either algorithm should
    # be the same
    result = PixelReplaceStep.call(input_model, skip=False, algorithm=algorithm)

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result.slits[0], ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result.slits[0], ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model.slits[0], ext)[bad_idx])
        assert getattr(result.slits[0], ext)[bad_idx] == 1.0

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.slits[0].dq[bad_idx] == (
        input_model.slits[0].dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.slits[0].dq[..., :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.slits[0].dq[..., 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    # Step is recorded as complete
    assert result.meta.cal_step.pixel_replace == "COMPLETE"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.pixel_replace is None

    result.close()
    input_model.close()


def test_pixel_replace_container_names(tmp_cwd, nirspec_fs_slitmodel):
    """Test pixel replace output names for input container."""
    input_model = nirspec_fs_slitmodel[0].copy()
    input_model.meta.filename = "jwst_nirspec_1_cal.fits"
    input_model2 = nirspec_fs_slitmodel[0].copy()
    input_model2.meta.filename = "jwst_nirspec_2_cal.fits"
    cfiles = [input_model, input_model2]
    container = ModelContainer(cfiles)

    expected_name = ["jwst_nirspec_1_pixelreplacestep.fits", "jwst_nirspec_2_pixelreplacestep.fits"]

    result = PixelReplaceStep.call(container, skip=False, save_results=True)
    for i, model in enumerate(result):
        assert model.meta.filename == expected_name[i]
        assert model.meta.cal_step.pixel_replace == "COMPLETE"

    result_files = glob(os.path.join(tmp_cwd, "*pixelreplacestep.fits"))
    for i, file in enumerate(sorted(result_files)):
        basename = os.path.basename(file)
        assert expected_name[i] == basename
        with datamodels.open(file) as model:
            assert model.meta.cal_step.pixel_replace == "COMPLETE"
            assert model.meta.filename == expected_name[i]

    # Input is not modified
    for model in container:
        assert model.meta.cal_step.pixel_replace is None

    result.close()
    input_model.close()


def test_pixel_replace_no_valid_data(caplog, nirspec_tso):
    """Test pixel replace for no valid data."""
    input_model, bad_idx = nirspec_tso

    # Set a middle region to NaN to test invalid data handling
    input_model.data[2, :, 5:15] = np.nan
    input_model.dq[2, 1:-1, 5:15] = 1

    # Set one pixel valid in the middle with no valid data next to
    # it to test missing adjacent data
    input_model.data[2, 10, 10] = 1.0
    input_model.dq[2, 10, 10] = 0
    input_model.dq[2, :, 9] = 1
    input_model.dq[2, :, 11] = 1

    result = PixelReplaceStep.call(input_model, algorithm="fit_profile", n_adjacent_cols=1)

    assert caplog.text.count("has no valid values - skipping") == 7
    assert caplog.text.count("has no valid adjacent values - skipping") == 1

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # non-science edges are uncorrected
        assert np.all(np.isnan(getattr(result, ext)[..., :, 1]))
        assert np.all(np.isnan(getattr(result, ext)[..., 1, :]))

        # bad pixel is replaced: input had one nan value, output does not
        assert np.isnan(getattr(input_model, ext)[bad_idx])
        assert getattr(result, ext)[bad_idx] == 1.0

        # invalid data is left alone
        assert np.all(np.isnan(result.data[2, 1:-1, 5:10]))
        assert np.all(np.isnan(result.data[2, 1:-1, 11:15]))

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated. The non-science edges are unchanged.
    assert result.dq[bad_idx] == (
        input_model.dq[bad_idx] - flags["DO_NOT_USE"] + flags["FLUX_ESTIMATED"]
    )
    assert np.all(result.dq[:2, :, 1] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])
    assert np.all(result.dq[:2, 1, :] == flags["DO_NOT_USE"] + flags["NON_SCIENCE"])

    # Invalid region is still marked DNU
    assert np.all(result.dq[2, 1:-1, 5:10] == flags["DO_NOT_USE"])
    assert np.all(result.dq[2, 1:-1, 11:15] == flags["DO_NOT_USE"])
    assert np.all(result.dq[2, 10, 10] == 0)
    assert np.all(result.dq[2, :, 9] == flags["DO_NOT_USE"])
    assert np.all(result.dq[2, :, 11] == flags["DO_NOT_USE"])

    result.close()
    input_model.close()


@pytest.mark.parametrize(
    "dataset", ["nirspec_tso", "miri_lrs", "miri_mrs", "nirspec_msa_multislit"]
)
@pytest.mark.parametrize("trace_present", [True, False, None])
def test_pixel_replace_with_trace_model(request, dataset, trace_present):
    """Test pixel replace with the trace model algorithm."""
    input_model, bad_idx = request.getfixturevalue(dataset)
    input_model = input_model.copy()

    # For this test, add a trace model with a different value at the bad index:
    # it should be ignored if trace_model is False.
    # Also set var_flat to None to make sure missing variance does not error.
    input_model.meta.cal_step.adaptive_trace_model = "COMPLETE"
    if "multislit" in dataset:
        slit = input_model.slits[0]
        if trace_present is not None:
            slit.trace_model = np.full(slit.data.shape[-2:], np.nan)
            if trace_present:
                slit.trace_model[bad_idx[-2:]] = 2.0
        slit.var_flat = None
    else:
        if trace_present is not None:
            input_model.trace_model = np.full(input_model.data.shape[-2:], np.nan)
            if trace_present:
                input_model.trace_model[bad_idx[-2:]] = 2.0
        input_model.var_flat = None

    result = PixelReplaceStep.call(input_model, algorithm="trace_model")

    if "multislit" in dataset:
        result = result.slits[0]

    for ext in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]:
        # bad pixel is replaced
        if ext == "var_flat":
            assert getattr(result, ext) is None
        elif ext == "data" and trace_present:
            # fixed via trace
            assert getattr(result, ext)[bad_idx] == 2.0
        else:
            # fixed via mingrad
            assert getattr(result, ext)[bad_idx] == 1.0

    # The DQ plane for the bad pixel is updated to remove do-not-use
    # and add flux-estimated.
    assert result.dq[bad_idx] == flags["OTHER_BAD_PIXEL"] + flags["FLUX_ESTIMATED"]


def test_pixel_replace_run_atm(nirspec_ifu):
    # only the IFU fixture has a proper WCS, so use that
    input_model, bad_idx = nirspec_ifu
    result = PixelReplaceStep.call(input_model, algorithm="trace_model")
    assert result.data[bad_idx] == 1.0
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"
    assert result.trace_model is not None


def test_pixel_replace_atm_previously_failed(caplog, nirspec_tso):
    input_model, bad_idx = nirspec_tso
    input_model = input_model.copy()
    input_model.meta.cal_step.adaptive_trace_model = "FAILED"
    result = PixelReplaceStep.call(input_model, algorithm="trace_model")

    # ATM is not called
    assert result.meta.cal_step.adaptive_trace_model == "FAILED"
    assert result.trace_model is None

    # Defaulting to mingrad fixes the bad pixel anyway
    assert "Defaulting to the 'mingrad' method" in caplog.text
    assert result.data[bad_idx] == 1.0


def test_pixel_replace_atm_error(caplog, nirspec_tso):
    input_model, bad_idx = nirspec_tso
    result = PixelReplaceStep.call(input_model, algorithm="trace_model")

    # ATM did not run
    assert result.meta.cal_step.adaptive_trace_model is None
    assert result.trace_model is None
    assert "ValueError: Unknown detector" in caplog.text

    # Defaulting to mingrad fixes the bad pixel anyway
    assert "Defaulting to the 'mingrad' method" in caplog.text
    assert result.data[bad_idx] == 1.0


def test_skip_unexpected_type():
    bad_model = datamodels.RampModel()
    result = PixelReplaceStep.call(bad_model)

    # Step is failed
    assert result.meta.cal_step.pixel_replace == "FAILED"

    # Input is not modified
    assert result is not bad_model
    assert bad_model.meta.cal_step.pixel_replace is None


def test_skip_unexpected_type_in_container():
    bad_model = datamodels.RampModel()
    container = ModelContainer([bad_model])
    result = PixelReplaceStep.call(container)

    # Step is failed
    assert result[0].meta.cal_step.pixel_replace == "FAILED"

    # Input is not modified
    assert result[0] is not bad_model
    assert bad_model.meta.cal_step.pixel_replace is None


@pytest.mark.parametrize("dataset", ["nirspec_tso", "nirspec_fs_slitmodel", "miri_mrs"])
@pytest.mark.parametrize("variance", ["var_poisson", "var_rnoise", "var_flat"])
@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_unset_variances(request, dataset, algorithm, variance):
    """
    Test that the step handles unset variance arrays gracefully.

    They should be left as None and raise no errors.
    """
    input_model = request.getfixturevalue(dataset)[0].copy()

    # Set variance to None
    input_model[variance] = None

    result = PixelReplaceStep.call(input_model, algorithm=algorithm)
    assert result[variance] is None


@pytest.mark.parametrize("variance", ["var_poisson", "var_rnoise", "var_flat"])
@pytest.mark.parametrize("algorithm", ["fit_profile", "mingrad"])
def test_unset_variances_multislit(nirspec_msa_multislit, algorithm, variance):
    """
    Test that the step handles unset variance arrays gracefully for multislit data.

    They should be left as None and raise no errors.
    """
    input_model = nirspec_msa_multislit[0].copy()

    # Set variance to None
    for slit in input_model.slits:
        setattr(slit, variance, None)

    result = PixelReplaceStep.call(input_model, algorithm=algorithm)

    # Variance is still None
    for slit in result.slits:
        assert getattr(slit, variance) is None


def test_n_replaced(caplog, nirspec_fs_slitmodel):
    input_model = nirspec_fs_slitmodel[0].copy()
    bad_idx = nirspec_fs_slitmodel[1]

    # Default call: DQ is DNU + other flag only
    assert input_model.dq[bad_idx] == flags["DO_NOT_USE"] | flags["OTHER_BAD_PIXEL"]
    result = PixelReplaceStep.call(input_model, algorithm="mingrad")
    assert "1 pixels replaced" in caplog.text
    assert result.dq[bad_idx] == flags["OTHER_BAD_PIXEL"] | flags["FLUX_ESTIMATED"]

    # Add flux estimated to the input flag: it should still be replaced and reported
    input_model.dq[bad_idx] |= flags["FLUX_ESTIMATED"]
    result = PixelReplaceStep.call(input_model, algorithm="mingrad")
    assert "1 pixels replaced" in caplog.text
    assert result.dq[bad_idx] == flags["OTHER_BAD_PIXEL"] | flags["FLUX_ESTIMATED"]

    # Calling again should report 0 replaced, since the bad pixel now has a finite value
    result = PixelReplaceStep.call(result, algorithm="mingrad")
    assert "0 pixels replaced" in caplog.text
    assert result.dq[bad_idx] == flags["OTHER_BAD_PIXEL"] | flags["FLUX_ESTIMATED"]


def test_missing_regions(nirspec_ifu):
    """
    Test fit_profile with an IFU file missing regions.

    This shouldn't happen for any current data, but could be encountered for old
    or custom-reduced data.
    """
    input_model = nirspec_ifu[0].copy()
    input_model.regions = None
    with pytest.raises(ValueError, match="missing region map"):
        PixelReplaceStep.call(input_model, algorithm="fit_profile")


def test_bad_algorithm(caplog, nirspec_fs_slitmodel):
    input_model = nirspec_fs_slitmodel[0]
    with pytest.raises(KeyError):
        PixelReplacement(input_model, algorithm="bad")
    assert "name 'bad' provided does not match" in caplog.text


def test_bad_input_model(caplog):
    input_model = datamodels.RampModel()
    pr = PixelReplacement(input_model)
    pr.replace()
    assert "not supported" in caplog.text


def test_no_good_pixels(caplog, nirspec_fs_slitmodel):
    input_model = nirspec_fs_slitmodel[0].copy()
    pr = PixelReplacement(input_model)
    arrays = pr._arrays_from_model(input_model)

    # all pixels bad
    arrays.dq |= flags["DO_NOT_USE"]

    result = pr.fit_profile(arrays)
    assert "No good pixels" in caplog.text
    assert np.all(result.dq & flags["DO_NOT_USE"] > 0)


def test_fit_profile_norm_scale(nirspec_fs_slitmodel):
    input_model = nirspec_fs_slitmodel[0].copy()
    bad_idx = nirspec_fs_slitmodel[1]
    pr = PixelReplacement(input_model)
    arrays = pr._arrays_from_model(input_model)

    # add non-unity values for flux and error to trigger profile scaling
    arrays.data *= 10
    arrays.err *= 0.1
    assert np.isnan(arrays.data[bad_idx])
    assert np.isnan(arrays.err[bad_idx])

    # profile scales correctly to fill in the expected value
    result = pr.fit_profile(arrays)
    assert np.isclose(result.data[bad_idx], 10)
    assert np.isclose(result.err[bad_idx], 0.1)
    assert result.dq[bad_idx] == flags["FLUX_ESTIMATED"] | flags["OTHER_BAD_PIXEL"]


def test_custom_slice():
    test_array = np.arange(100).reshape((10, 10))
    pr = PixelReplacement(None)

    # middle x
    horiz_slice = pr.custom_slice(1, 5)
    np.testing.assert_equal(test_array[horiz_slice], np.arange(5, 100, 10))

    # middle y
    vert_slice = pr.custom_slice(2, 5)
    np.testing.assert_equal(test_array[vert_slice], np.arange(50, 60, 1))

    # invalid
    with pytest.raises(IndexError, match="requires valid"):
        pr.custom_slice(3, 5)


@pytest.mark.parametrize("dispersion_direction", [1, 2])
@pytest.mark.parametrize("bad_value", [0, np.nan])
@pytest.mark.parametrize(
    "xindx,yindx",
    [
        (1, 2),
        ([3, 5], [6, 7]),
        ([3, 3, 3], [4, 5, 6]),
        ([4, 5, 6], [3, 3, 3]),
        ([4, 4, 5, 5], [4, 5, 4, 5]),
    ],
)
def test_interp_along_wavelength(dispersion_direction, bad_value, xindx, yindx):
    pr = PixelReplacement(None)

    data = np.tile(np.arange(10, dtype=float), 10).reshape(10, 10)
    expected = data.copy()

    xindx = np.array(xindx)
    yindx = np.array(yindx)
    data[yindx, xindx] = bad_value
    pr._interp_along_wavelength(data, dispersion_direction, xindx, yindx)
    np.testing.assert_allclose(data, expected)


def test_interp_along_wavelength_no_valid():
    data = np.full((10, 10), np.nan)
    xindx = np.array([3, 5])
    yindx = np.array([6, 7])
    pr = PixelReplacement(None)
    pr._interp_along_wavelength(data, 1, xindx, yindx)
    assert np.all(np.isnan(data))


def test_interp_along_wavelength_edge_value():
    data = np.tile(np.arange(10, dtype=float), 10).reshape(10, 10)
    expected = data.copy()

    xindx = np.array([9])
    yindx = np.array([9])
    data[9, 9] = np.nan

    # Edge pixel is extended from last available
    expected[9, 9] = 8

    pr = PixelReplacement(None)
    pr._interp_along_wavelength(data, 1, xindx, yindx)
    np.testing.assert_allclose(data, expected)
