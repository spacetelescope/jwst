import json

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.adaptive_trace_model.adaptive_trace_model_step import AdaptiveTraceModelStep
from jwst.adaptive_trace_model.tests import helpers
from jwst.datamodels import ModelContainer


@pytest.fixture(scope="module")
def miri_mrs_model():
    model = helpers.miri_mrs_model()
    yield model
    model.close()


@pytest.fixture(scope="module")
def nirspec_ifu_model_with_source():
    model = helpers.nirspec_ifu_model_with_source()
    yield model
    model.close()


@pytest.fixture(scope="module")
def miri_mrs_model_with_source():
    model = helpers.miri_mrs_model_with_source()
    yield model
    model.close()


@pytest.fixture(scope="module")
def nirspec_ifu_slice_wcs():
    model = helpers.nirspec_ifu_model_with_source(wcs_style="slice")
    yield model
    model.close()


@pytest.fixture()
def asn_input(tmp_path, miri_mrs_model):
    """
    Create an association file with two MIRI MRS inputs.

    The science images are saved to a temporary directory so
    the association can locate them when loaded.

    Returns
    -------
    asn : dict
        The association dictionary.
    """
    models = [miri_mrs_model.copy(), miri_mrs_model.copy()]
    filenames = ["test1_cal.fits", "test2_cal.fits"]
    for model, filename in zip(models, filenames):
        model.save(str(tmp_path / filename))
        model.close()

    asn = {
        "asn_type": "test",
        "asn_pool": "test",
        "asn_id": "o001",
        "products": [
            {
                "name": "product_a",
                "members": [
                    {"expname": filenames[0], "exptype": "science"},
                    {"expname": filenames[1], "exptype": "science"},
                ],
            },
        ],
    }

    return asn


@pytest.mark.parametrize("dataset", ["miri_mrs_model", "nirspec_ifu_model_with_source"])
def test_adaptive_trace_model_step_success(request, dataset):
    model = request.getfixturevalue(dataset)

    # run with a high threshold so spline fits are not performed, for speed
    result = AdaptiveTraceModelStep.call(model, oversample=1, fit_threshold=100000)

    # step is complete
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"

    # output is not the same as input
    assert result is not model
    assert model.meta.cal_step.adaptive_trace_model is None

    # data is unchanged with oversample=1
    np.testing.assert_equal(result.data, model.data)

    # trace is attached, but contains NaN only for high threshold
    assert result.trace_model.shape == result.data.shape
    assert np.all(np.isnan(result.trace_model))

    result.close()


def test_adaptive_trace_model_step_with_source(miri_mrs_model_with_source):
    model = miri_mrs_model_with_source
    result = AdaptiveTraceModelStep.call(model, oversample=1, slope_limit=0)
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"

    # data is unchanged with oversample=1
    np.testing.assert_equal(result.data, model.data)

    # trace is attached, contains non-NaN trace for the one bright slit only
    assert result.trace_model.shape == result.data.shape
    det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
    indx = det2ab_transform.label_mapper.mapper == 120
    assert np.all(np.isnan(result.trace_model[~indx]))
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.85 * np.sum(indx)

    # fit trace is a reasonable model of the slice but not perfect
    valid = indx & ~np.isnan(result.data) & ~np.isnan(result.trace_model)
    atol = 0.25 * np.nanmax(model.data)
    np.testing.assert_allclose(result.data[valid], result.trace_model[valid], atol=atol)
    result.close()


def test_adaptive_trace_model_step_negative_mean(miri_mrs_model_with_source):
    model = miri_mrs_model_with_source.copy()

    # Add a significant negative mean value to the data
    original_max = np.nanmax(model.data)
    model.data -= 100

    result = AdaptiveTraceModelStep.call(model, oversample=1, slope_limit=0)
    det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
    indx = det2ab_transform.label_mapper.mapper == 120
    assert np.all(np.isnan(result.trace_model[~indx]))
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.85 * np.sum(indx)

    # fit trace is a reasonable model of the slice, minus the negative mean
    valid = indx & ~np.isnan(result.data) & ~np.isnan(result.trace_model)
    atol = 0.25 * original_max
    np.testing.assert_allclose(result.data[valid] + 100, result.trace_model[valid], atol=atol)

    model.close()
    result.close()


def test_adaptive_trace_model_step_oversample(miri_mrs_model):
    model = miri_mrs_model
    result = AdaptiveTraceModelStep.call(model, oversample=2)
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"

    # data is twice the size of the input along the x axis
    extnames = ["data", "dq", "err", "var_poisson", "var_rnoise", "var_flat"]
    for extname in extnames:
        input_ext = getattr(model, extname)
        output_ext = getattr(result, extname)
        assert output_ext.shape == (input_ext.shape[0], 2 * input_ext.shape[1])

        # Check mean value: errors are inflated, data should match
        inflation = 0.23 * 2 + 0.77
        if extname == "data":
            np.testing.assert_allclose(np.nanmean(output_ext), np.mean(input_ext), atol=1e-7)
        elif extname == "err":
            np.testing.assert_allclose(
                np.nanmean(output_ext), np.mean(input_ext) * inflation, atol=1e-7
            )
        elif extname.startswith("var"):
            np.testing.assert_allclose(
                np.nanmean(output_ext), np.mean(input_ext) * inflation**2, atol=1e-7
            )

    # trace is attached, but contains NaN only for flat input data
    assert result.trace_model.shape == result.data.shape
    assert np.all(np.isnan(result.trace_model))

    result.close()


@pytest.mark.parametrize(
    "dataset",
    ["miri_mrs_model_with_source", "nirspec_ifu_model_with_source", "nirspec_ifu_slice_wcs"],
)
def test_adaptive_trace_model_step_oversample_with_source(request, dataset):
    model = request.getfixturevalue(dataset)

    fit_threshold = 10.0
    slope_limit = 0.05
    result = AdaptiveTraceModelStep.call(
        model, oversample=2, slope_limit=slope_limit, fit_threshold=fit_threshold
    )
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"

    # data is twice the size of the input along the x axis
    extnames = ["data", "dq", "err", "var_poisson", "var_rnoise", "var_flat"]
    for extname in extnames:
        # check for extension presence
        if not model.hasattr(extname):
            assert not result.hasattr(extname)
            continue

        input_ext = getattr(model, extname)
        output_ext = getattr(result, extname)
        if dataset.startswith("mir"):
            assert output_ext.shape == (input_ext.shape[0], 2 * input_ext.shape[1])
        else:
            assert output_ext.shape == (input_ext.shape[0] * 2, input_ext.shape[1])

    # trace is attached, contains non-NaN trace for the one bright slit only
    assert result.trace_model.shape == result.data.shape
    if dataset.startswith("mir"):
        indx = result.regions == 120
    else:
        indx = result.regions == 16
    assert np.all(np.isnan(result.trace_model[~indx]))
    # Only the compact core is modeled, so ~20% of the slice is non-NaN
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.2 * np.sum(indx)

    # fit trace is a reasonable model of the slice but not perfect -
    # the slice is mostly linearly interpolated
    valid = indx & ~np.isnan(result.data) & ~np.isnan(result.trace_model)
    atol = 0.25 * np.nanmax(model.data)
    np.testing.assert_allclose(result.data[valid], result.trace_model[valid], atol=atol)

    result.close()


def test_adaptive_trace_model_step_with_container(miri_mrs_model):
    model = miri_mrs_model
    container = ModelContainer([model, model.copy()])
    result = AdaptiveTraceModelStep.call(container, oversample=1)

    assert result is not container
    assert isinstance(result, ModelContainer)
    for input_model, output_model in zip(container, result, strict=True):
        assert input_model is not output_model
        assert input_model.meta.cal_step.adaptive_trace_model is None
        assert output_model.meta.cal_step.adaptive_trace_model == "COMPLETE"
        assert output_model.hasattr("trace_model")


def test_adaptive_trace_model_unsupported_model(caplog):
    model = datamodels.ImageModel()
    result = AdaptiveTraceModelStep.call(model)
    assert "only implemented for IFU" in caplog.text

    assert result is not model
    assert result.meta.cal_step.adaptive_trace_model == "SKIPPED"
    assert model.meta.cal_step.adaptive_trace_model is None


@pytest.mark.slow
def test_adaptive_trace_model_step_psf_optimal(caplog, miri_mrs_model_with_source):
    model = miri_mrs_model_with_source
    result = AdaptiveTraceModelStep.call(model, oversample=2, psf_optimal=True)
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"
    assert "Ignoring fit threshold and slope limit" in caplog.text

    # trace is attached, contains non-NaN trace for all slices
    assert result.trace_model.shape == result.data.shape
    indx = result.regions > 0
    assert np.all(np.isnan(result.trace_model[~indx]))
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.9 * np.sum(indx)

    # fit trace is exactly the same as the data, since the spline
    # model is used everywhere
    valid = indx & ~np.isnan(result.data)
    np.testing.assert_allclose(result.data[valid], result.trace_model[valid])

    result.close()


def test_save_container_asn_id_present(tmp_path, asn_input):
    asn_file = str(tmp_path / "test_asn.json")
    with open(asn_file, "w") as fh:
        json.dump(asn_input, fh)

    AdaptiveTraceModelStep.call(asn_file, output_dir=str(tmp_path), suffix="atm", save_results=True)
    expected_output = [tmp_path / "test1_o001_atm.fits", tmp_path / "test2_o001_atm.fits"]
    for output in expected_output:
        assert output.exists()


def test_save_container_asn_id_missing(tmp_path, asn_input):
    del asn_input["asn_id"]
    asn_file = str(tmp_path / "test_asn.json")
    with open(asn_file, "w") as fh:
        json.dump(asn_input, fh)

    AdaptiveTraceModelStep.call(asn_file, output_dir=str(tmp_path), suffix="atm", save_results=True)
    expected_output = [tmp_path / "test1_atm.fits", tmp_path / "test2_atm.fits"]
    for output in expected_output:
        assert output.exists()


@pytest.mark.parametrize("oversample", [1.0, 2.0])
def test_adaptive_trace_model_step_save_intermediate(tmp_path, miri_mrs_model, oversample):
    model = miri_mrs_model
    AdaptiveTraceModelStep.call(
        model,
        oversample=oversample,
        save_results=True,
        suffix="atm",
        save_intermediate_results=True,
        output_dir=str(tmp_path),
    )

    # Check for expected files
    expected = [
        "test12SHORT_atm.fits",
        "test12SHORT_spline_full.fits",
        "test12SHORT_spline_used.fits",
    ]
    expect_empty = [False, True, True]
    if oversample > 1:
        # Extra files expected if oversampling is done
        expected.extend(["test12SHORT_linear_interp.fits", "test12SHORT_spline_residual.fits"])
        expect_empty.extend([False, True])
    for filename, is_empty in zip(expected, expect_empty):
        assert (tmp_path / filename).exists()
        with datamodels.open(str(tmp_path / filename)) as model:
            assert isinstance(model, datamodels.IFUImageModel)
            if is_empty:
                assert np.all(np.isnan(model.data))
            else:
                assert not np.all(np.isnan(model.data))
