import numpy as np
import pytest
from stdatamodels.jwst.datamodels import ImageModel

from jwst.adaptive_trace_model.adaptive_trace_model_step import AdaptiveTraceModelStep
from jwst.adaptive_trace_model.tests import helpers
from jwst.datamodels import ModelContainer


@pytest.mark.parametrize("mode", ["MIR_MRS", "NRS_IFU"])
def test_adaptive_trace_model_step_success(mode):
    if mode == "MIR_MRS":
        model = helpers.miri_mrs_model()
    else:
        model = helpers.nirspec_ifu_model_with_source()
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

    model.close()
    result.close()


def test_adaptive_trace_model_step_with_source():
    model = helpers.miri_mrs_model_with_source()
    result = AdaptiveTraceModelStep.call(model, oversample=1)
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"

    # data is unchanged with oversample=1
    np.testing.assert_equal(result.data, model.data)

    # trace is attached, contains non-NaN trace for the one bright slit only
    assert result.trace_model.shape == result.data.shape
    det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
    indx = det2ab_transform.label_mapper.mapper == 120
    assert np.all(np.isnan(result.trace_model[~indx]))
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.9 * np.sum(indx)

    # fit trace is a reasonable model of the slice but not perfect
    valid = indx & ~np.isnan(result.data)
    atol = 0.25 * np.nanmax(model.data)
    np.testing.assert_allclose(result.data[valid], result.trace_model[valid], atol=atol)

    model.close()
    result.close()


def test_adaptive_trace_model_step_negative_mean():
    model = helpers.miri_mrs_model_with_source()

    # Add a significant negative mean value to the data
    original_max = np.nanmax(model.data)
    model.data -= 100

    result = AdaptiveTraceModelStep.call(model, oversample=1)
    det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
    indx = det2ab_transform.label_mapper.mapper == 120
    assert np.all(np.isnan(result.trace_model[~indx]))
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.9 * np.sum(indx)

    # fit trace is a reasonable model of the slice, minus the negative mean
    valid = indx & ~np.isnan(result.data)
    atol = 0.25 * original_max
    np.testing.assert_allclose(result.data[valid] + 100, result.trace_model[valid], atol=atol)

    model.close()
    result.close()


def test_adaptive_trace_model_step_oversample():
    model = helpers.miri_mrs_model()
    result = AdaptiveTraceModelStep.call(model, oversample=2)
    assert result.meta.cal_step.adaptive_trace_model == "COMPLETE"

    # data is twice the size of the input along the x axis
    extnames = ["data", "dq", "err", "var_poisson", "var_rnoise", "var_flat"]
    for extname in extnames:
        input_ext = getattr(model, extname)
        output_ext = getattr(result, extname)
        assert output_ext.shape == (input_ext.shape[0], 2 * input_ext.shape[1])
        if extname != "dq":
            np.testing.assert_allclose(np.nanmean(output_ext), np.mean(input_ext), atol=1e-7)

    # trace is attached, but contains NaN only for flat input data
    assert result.trace_model.shape == result.data.shape
    assert np.all(np.isnan(result.trace_model))

    model.close()
    result.close()


@pytest.mark.parametrize("mode", ["MIR_MRS", "NRS_IFU", "NRS_IFU_SLICE_WCS"])
def test_adaptive_trace_model_step_oversample_with_source(mode):
    fit_threshold = 10.0
    if mode == "MIR_MRS":
        model = helpers.miri_mrs_model_with_source()
    elif mode == "NRS_IFU_SLICE_WCS":
        model = helpers.nirspec_ifu_model_with_source(wcs_style="slice")
    else:
        model = helpers.nirspec_ifu_model_with_source()
    result = AdaptiveTraceModelStep.call(
        model, oversample=2, slope_limit=0.05, fit_threshold=fit_threshold
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
        if mode.startswith("MIR"):
            assert output_ext.shape == (input_ext.shape[0], 2 * input_ext.shape[1])
        else:
            assert output_ext.shape == (input_ext.shape[0] * 2, input_ext.shape[1])

    # trace is attached, contains non-NaN trace for the one bright slit only
    assert result.trace_model.shape == result.data.shape
    if mode == "MIR_MRS":
        indx = result.regions == 120
    else:
        indx = result.regions == 16
    assert np.all(np.isnan(result.trace_model[~indx]))
    assert np.sum(~np.isnan(result.trace_model[indx])) > 0.9 * np.sum(indx)

    # fit trace is a reasonable model of the slice but not perfect -
    # the slice is mostly linearly interpolated
    valid = indx & ~np.isnan(result.data)
    atol = 0.25 * np.nanmax(model.data)
    np.testing.assert_allclose(result.data[valid], result.trace_model[valid], atol=atol)

    model.close()
    result.close()


def test_adaptive_trace_model_step_with_container():
    model = helpers.miri_mrs_model()
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
    model = ImageModel()
    result = AdaptiveTraceModelStep.call(model)
    assert "only implemented for IFU" in caplog.text

    assert result is not model
    assert result.meta.cal_step.adaptive_trace_model == "SKIPPED"
    assert model.meta.cal_step.adaptive_trace_model is None
