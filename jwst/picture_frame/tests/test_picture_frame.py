import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise.tests.helpers import make_nirspec_fs_model, make_nrs_fs_full_ramp
from jwst.lib.basic_utils import LoggingContext
from jwst.picture_frame import picture_frame as pf
from jwst.picture_frame.tests.helpers import picture_frame_model


def test_correction_rate(caplog):
    pctfrm_model = picture_frame_model()
    input_model = make_nirspec_fs_model()

    # Add a edge region to correct
    input_model.data[
        pf.CENTER_REGION[0] : pf.CENTER_REGION[1], pf.CENTER_REGION[0] : pf.CENTER_REGION[1]
    ] += 0.2

    # Input data is modified in place, so work on a copy
    with LoggingContext(logging.getLogger("jwst"), level=logging.DEBUG):
        cleaned, mask, correction, status = pf.correct_picture_frame(
            input_model.copy(), pctfrm_model, save_correction=True, save_mask=True
        )

    assert "Median center reference: 1.1" in caplog.text
    assert "Median edge reference: 0.1" in caplog.text
    assert "Median center: 1.2" in caplog.text
    assert "Median edge: 1.0" in caplog.text

    # Flat artifact should be perfectly corrected
    np.testing.assert_allclose(cleaned.data, 0.0)

    # Correction data should match input type and contain output - input data
    assert isinstance(mask, datamodels.ImageModel)
    np.testing.assert_allclose(correction.data, cleaned.data - input_model.data)

    # Mask should not be flat, since science regions are masked
    assert isinstance(mask, datamodels.ImageModel)
    assert not np.all(mask.data == 1)


def test_correction_rateints():
    pctfrm_model = picture_frame_model()

    # Make a rateints model by repeating a rate image a couple times
    rate_model = make_nirspec_fs_model()
    input_model = datamodels.CubeModel((3, *rate_model.data.shape))
    input_model.update(rate_model)
    for i in range(3):
        input_model.data[i] = rate_model.data + i

    # Add a edge region to correct
    input_model.data[
        :, pf.CENTER_REGION[0] : pf.CENTER_REGION[1], pf.CENTER_REGION[0] : pf.CENTER_REGION[1]
    ] += 0.2

    # Input data is modified in place, so work on a copy
    cleaned, mask, correction, status = pf.correct_picture_frame(
        input_model.copy(), pctfrm_model, save_correction=True, save_mask=True
    )

    # Flat artifact should be perfectly corrected
    np.testing.assert_allclose(cleaned.data, 0.0)

    # Correction data should match input type and contain output - input data
    assert isinstance(correction, datamodels.CubeModel)
    np.testing.assert_allclose(correction.data, cleaned.data - input_model.data)

    # Mask should not be flat, since science regions are masked
    assert isinstance(mask, datamodels.ImageModel)
    assert not np.all(mask.data == 1)


def test_correction_ramp():
    pctfrm_model = picture_frame_model()
    input_model = make_nrs_fs_full_ramp()

    # Add a edge region
    input_model.data[
        :, :, pf.CENTER_REGION[0] : pf.CENTER_REGION[1], pf.CENTER_REGION[0] : pf.CENTER_REGION[1]
    ] += 0.2

    # Input data is modified in place, so work on a copy
    cleaned, mask, correction, status = pf.correct_picture_frame(
        input_model.copy(), pctfrm_model, save_correction=True, save_mask=True
    )

    # The first group is unchanged, the second one is cleaned to match the first
    np.testing.assert_allclose(cleaned.data[0, 0], input_model.data[0, 0])
    np.testing.assert_allclose(cleaned.data[0, 1], input_model.data[0, 0])

    # Correction data should match input type and contain output - input data
    assert isinstance(correction, datamodels.RampModel)
    np.testing.assert_allclose(correction.data, cleaned.data - input_model.data)

    # Mask should not be flat, since science regions are masked
    assert isinstance(mask, datamodels.ImageModel)
    assert not np.all(mask.data == 1)


def test_correction_invalid_subarray(caplog):
    pctfrm_model = picture_frame_model()

    # Modify data to have an invalid subarray
    input_model = make_nirspec_fs_model()
    input_model.meta.subarray.name = "ALLSLITS"

    cleaned, mask, correction, status = pf.correct_picture_frame(
        input_model.copy(), pctfrm_model, save_correction=True, save_mask=True
    )

    # Message is logged, mask and correction are None
    assert "only applicable to NIRSpec full frame" in caplog.text
    np.testing.assert_allclose(cleaned.data, input_model.data)
    assert status == "SKIPPED"
    assert mask is None
    assert correction is None


def test_correction_invalid_data(caplog):
    pctfrm_model = picture_frame_model()

    # Modify data to have all invalid values
    input_model = make_nirspec_fs_model()
    input_model.data *= np.nan

    cleaned, mask, correction, status = pf.correct_picture_frame(
        input_model.copy(), pctfrm_model, save_correction=True, save_mask=True
    )

    # Message is logged, correction is None
    assert "No data to scale" in caplog.text
    np.testing.assert_allclose(cleaned.data, input_model.data)
    assert status == "SKIPPED"
    assert correction is None

    # Mask is calculated and returned: it shows no valid background data to fit
    assert mask is not None
    assert np.all(mask.data == 0)
