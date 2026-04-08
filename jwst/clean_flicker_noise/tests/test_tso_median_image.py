import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.clean_flicker_noise import tso_median_image as tmi
from jwst.clean_flicker_noise.background_level import background_level
from jwst.clean_flicker_noise.tests import helpers


def test_make_median_ndim_2():
    input_model = helpers.make_small_rate_model()
    with pytest.raises(ValueError, match="2D data"):
        tmi.make_median_image(input_model, input_model)

    input_model.close()


def test_make_median_rateints_one_int():
    input_model = helpers.make_small_rateints_model()
    shape = input_model.data.shape
    input_model.data = input_model.data[0].reshape((1, *shape[1:]))
    with pytest.raises(ValueError, match="<2 integrations"):
        tmi.make_median_image(input_model, input_model)

    input_model.close()


def test_make_median_rateints_non_tso():
    # Test a non-TSO mode
    input_model = helpers.make_small_rateints_model()
    med_img = tmi.make_median_image(input_model, input_model)
    assert med_img.shape == input_model.data.shape

    # Flat data, so the median image should be the same as the input data
    np.testing.assert_allclose(med_img, input_model.data)

    input_model.close()


def test_make_median_ramp_non_tso():
    # Test a non-TSO mode
    input_model = helpers.make_small_ramp_model()
    rateints_model = helpers.make_small_rateints_model()
    med_img = tmi.make_median_image(input_model, rateints_model)
    assert med_img.shape == input_model.data.shape

    # Flat data, so the median image should be the same as the input data
    np.testing.assert_allclose(med_img, input_model.data)

    input_model.close()


def test_non_tso_all_invalid():
    input_model = helpers.make_small_rateints_model()
    input_model.data *= np.nan
    with pytest.raises(ValueError, match="No valid flux for scaling"):
        tmi.make_median_image(input_model, input_model)

    input_model.close()


def test_non_tso_some_invalid():
    input_model = helpers.make_small_rateints_model()
    input_model.data[1] *= np.nan
    med_img = tmi.make_median_image(input_model, input_model)

    # invalid values replaced with median in output image
    np.testing.assert_allclose(med_img, 1.0)

    input_model.close()


def test_soss_extract():
    # Input data is filled with ones
    input_model = helpers.make_niriss_soss_rateints()
    nint = input_model.data.shape[0]

    # Extract a simple summed spectrum from order 0 with box width 15
    tso_spec = tmi._soss_box_extract(input_model)
    assert isinstance(tso_spec, datamodels.TSOMultiSpecModel)
    assert len(tso_spec.spec) == 1

    # Output spectrum has wavelength, summed flux, mocked integration time
    table = tso_spec.spec[0].spec_table
    np.testing.assert_allclose(table["FLUX"], 15)
    np.testing.assert_allclose(table["WAVELENGTH"].min(), 0.85, atol=0.1)
    np.testing.assert_allclose(table["WAVELENGTH"].max(), 2.83, atol=0.1)
    np.testing.assert_allclose(table["MJD-AVG"], np.arange(nint), atol=0.1)

    input_model.close()


def test_make_background_ramp():
    # Input ramp with non-trivial readout pattern
    ramp_shape = (2, 3, 4, 5)
    input_model = datamodels.RampModel(ramp_shape)
    input_model.meta.exposure.nframes = 2
    input_model.meta.exposure.groupgap = 3
    input_model.meta.exposure.frame_time = 4

    # Input background rateints with non-trivial data
    rng = np.random.default_rng(seed=42)
    bg_data = rng.normal(0, 1, (2, 4, 5))

    bg_ramp = tmi._make_background_ramp(input_model, bg_data)
    assert bg_ramp.shape == ramp_shape
    for i, integ in enumerate(bg_ramp):
        for j, group in enumerate(integ):
            # each group is rate * (nframes + (nframes + groupgap) * i) * frame_time
            assert np.allclose(group, bg_data[i] * (2 + 5 * j) * 4)

    input_model.close()


def test_soss_background_all_invalid():
    input_model = helpers.make_niriss_soss_rateints()
    input_model.data *= np.nan
    with pytest.raises(ValueError, match="No valid values in background rate"):
        tmi.make_median_image(input_model, input_model)
    input_model.close()


@pytest.mark.parametrize("model_failed", [True, False])
def test_soss_background_replace_invalid(monkeypatch, model_failed):
    input_model = helpers.make_niriss_soss_rateints()

    # Make a NaN value in each integration. Input values are otherwise all 1.0
    for i in range(input_model.data.shape[0]):
        input_model.data[i, 10 + i * 10, 10 + i * 10] = np.nan

    # Mock a failure in the model-based interpolation: it will use the median value instead
    if model_failed:

        def mock_bg_level(*args, **kwargs):
            return background_level(*args, background_method="median")

        monkeypatch.setattr(tmi, "background_level", mock_bg_level)

    # Invalid values are replaced with interpolated values in the output image
    med_image = tmi.make_median_image(input_model, input_model)
    np.testing.assert_allclose(med_image, 1.0)

    input_model.close()


def test_nrs_bots_extraction_fail():
    # Input data with no WCS: extraction will fail
    input_model = helpers.make_nrs_bots_rateints()
    with pytest.raises(ValueError, match="No valid spectra"):
        tmi.make_median_image(input_model, input_model)

    input_model.close()


def test_nrs_bots_success():
    # Smoke test NIRSpec BOTS with flat data
    input_model = helpers.make_nrs_bots_rateints()
    wcs_model = AssignWcsStep.call(input_model)

    med_image = tmi.make_median_image(wcs_model, wcs_model)
    assert med_image.shape == wcs_model.data.shape
    np.testing.assert_allclose(med_image, 1.0)

    input_model.close()
    wcs_model.close()


def test_miri_image_tso():
    # Smoke test MIRI imaging TSO with flat data
    input_model = helpers.make_miri_image_tso_rateints()
    wcs_model = AssignWcsStep.call(input_model)

    med_image = tmi.make_median_image(wcs_model, wcs_model)
    assert med_image.shape == wcs_model.data.shape
    np.testing.assert_allclose(med_image, 1.0)

    input_model.close()
    wcs_model.close()


@pytest.mark.parametrize("use_ramp", [True, False])
def test_input_not_modified(use_ramp):
    rateints_model = helpers.make_small_rateints_model()
    if use_ramp:
        input_model = helpers.make_small_ramp_model()
    else:
        input_model = rateints_model.copy()

    # add some noise to the input to make the median meaningful
    rng = np.random.default_rng(seed=42)
    input_model.data += rng.normal(0, 0.01, input_model.data.shape)
    rateints_model.data += rng.normal(0, 0.01, rateints_model.data.shape)

    input_copy = input_model.data.copy()
    rate_copy = rateints_model.data.copy()

    med_img = tmi.make_median_image(input_model, rateints_model)
    assert med_img.shape == input_model.data.shape

    # The median image should be close to the input data
    np.testing.assert_allclose(med_img, input_model.data, atol=0.05)

    # Input is not modified
    np.testing.assert_array_equal(input_model.data, input_copy)
    np.testing.assert_array_equal(rateints_model.data, rate_copy)

    input_model.close()
    rateints_model.close()
