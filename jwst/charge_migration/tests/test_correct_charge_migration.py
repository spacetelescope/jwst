import numpy as np
import pytest

from jwst.charge_migration.charge_migration import correct_charge_migration
from jwst.charge_migration.tests import helpers

IMG_SHAPE = (9, 9)
THRESHOLD = 100.0
PIXEL = (3, 4)


@pytest.fixture
def small_ramp():
    # Function scoped to allow update in place - model is small
    shape = (1, 3, *IMG_SHAPE)
    ramp_model = helpers.niriss_ramp_model_with_signal(THRESHOLD, PIXEL, shape=shape, start_group=1)
    return ramp_model


@pytest.fixture
def sigclip_ramp():
    # Function scoped to allow update in place - model is small
    shape = (2, 10, *IMG_SHAPE)
    ramp_model = helpers.niriss_ramp_model_with_signal(THRESHOLD, PIXEL, shape=shape, start_group=3)
    return ramp_model


@pytest.fixture(scope="module")
def readnoise_model():
    return helpers.mock_readnoise_model(IMG_SHAPE)


def test_correct_one_pix(caplog, small_ramp, readnoise_model):
    caplog.set_level("INFO", "jwst")

    # keep the input values at the expected pixel
    input_values = small_ramp.data.copy()

    # correct the ramp
    result = correct_charge_migration(
        small_ramp, readnoise_model, correction_threshold=THRESHOLD, mask_pad=0
    )

    # output is updated in place
    assert result is small_ramp

    # found 1 pixel to correct with mask_pad = 0
    assert "Identified 1 pixels to fit" in caplog.text

    # last group in expected pixel is corrected to a higher value (charge is restored)
    corrected = (0, 2, *PIXEL)
    assert result.data[corrected] > input_values[corrected]

    # other data is unaffected
    flat_index = np.ravel_multi_index(corrected, input_values.shape)
    np.testing.assert_allclose(
        np.delete(result.data, flat_index), np.delete(input_values, flat_index)
    )


def test_correct_neighbors(caplog, small_ramp, readnoise_model):
    caplog.set_level("INFO", "jwst")
    input_values = small_ramp.data.copy()

    # set fit quality limits to zero so all fits are used
    result, fit_info = correct_charge_migration(
        small_ramp,
        readnoise_model,
        correction_threshold=THRESHOLD,
        fit_info=True,
        snr_limit=0,
        rms_limit=0,
        mask_pad=2,
    )

    # 25 pixels are fit: expected one at center + all neighbors within mask_pad
    assert "Identified 25 pixels to fit" in caplog.text
    assert len(fit_info) == 25
    assert np.all(fit_info["fit_success"])

    # expected pixel is corrected, others unaffected: actual fit is zero for flat ramps
    corrected = (0, 2, *PIXEL)
    flat_index = np.ravel_multi_index(corrected, input_values.shape)
    assert result.data[corrected] > input_values[corrected]
    np.testing.assert_allclose(
        np.delete(result.data, flat_index), np.delete(input_values, flat_index)
    )


def test_clip_outliers(caplog, sigclip_ramp, readnoise_model):
    caplog.set_level("INFO", "jwst")
    input_values = sigclip_ramp.data.copy()
    result = correct_charge_migration(sigclip_ramp, readnoise_model, correction_threshold=THRESHOLD)

    # input data has one diff non-zero: it is sigma clipped out, so no correction is applied
    np.testing.assert_allclose(result.data, input_values)
    assert "No pixels to correct after sigma clipping" in caplog.text


def test_no_valid(caplog, small_ramp, readnoise_model):
    caplog.set_level("INFO", "jwst")
    small_ramp.data *= np.inf
    input_values = small_ramp.data.copy()
    result = correct_charge_migration(small_ramp, readnoise_model, correction_threshold=THRESHOLD)

    # input data has all invalid values: no correction is applied
    np.testing.assert_allclose(result.data, input_values)
    assert "No valid pixels to correct" in caplog.text


def test_no_pix_above_threshold(caplog, small_ramp, readnoise_model):
    caplog.set_level("INFO", "jwst")

    input_values = small_ramp.data.copy()
    result = correct_charge_migration(small_ramp, readnoise_model, correction_threshold=1e6)

    # input data is all below threshold: no correction is applied
    np.testing.assert_allclose(result.data, input_values)
    assert "No pixels over threshold to correct" in caplog.text


def test_too_few_diffs(caplog, small_ramp, readnoise_model):
    caplog.set_level("INFO", "jwst")

    # keep just 1 diff, so correction will bail
    small_ramp.data = small_ramp.data[:, :2]
    small_ramp.groupdq = small_ramp.groupdq[:, :2]

    input_values = small_ramp.data.copy()
    result = correct_charge_migration(small_ramp, readnoise_model, correction_threshold=THRESHOLD)

    # input data has too few diffs: no correction is applied
    np.testing.assert_allclose(result.data, input_values)
    assert "No valid pixels to correct" in caplog.text
