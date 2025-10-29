import os

import pytest
from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise import CleanFlickerNoiseStep, autoparam
from jwst.clean_flicker_noise.tests.helpers import (
    make_nircam_rate_model,
    make_niriss_rate_model,
    make_nirspec_fs_model,
    make_small_ramp_model,
)
from jwst.pipeline import Detector1Pipeline


@pytest.mark.parametrize("skip", [True, False])
def test_output_type(skip):
    input_model = make_small_ramp_model()
    step_instance = CleanFlickerNoiseStep()
    # Default is that step is skipped
    assert step_instance.skip == True
    cleaned = step_instance.call(input_model, skip=skip)

    # output is a ramp model either way
    assert isinstance(cleaned, datamodels.RampModel)

    input_model.close()
    cleaned.close()


@pytest.mark.parametrize("skip", [True, False])
def test_run_in_pipeline(skip):
    input_model = make_small_ramp_model()
    pipeline_instance = Detector1Pipeline()

    assert pipeline_instance.steps["clean_flicker_noise"]["skip"] == True

    # Run the pipeline, omitting steps that are incompatible with our datamodel
    cleaned = pipeline_instance.call(
        input_model,
        steps={
            "clean_flicker_noise": {"skip": skip},
            "ipc": {"skip": True},
            "reset": {"skip": True},
            "dark_current": {"skip": True},
        },
    )

    if skip:
        assert cleaned.meta.cal_step.clean_flicker_noise == "SKIPPED"
    else:
        assert cleaned.meta.cal_step.clean_flicker_noise == "COMPLETE"

    # Either way, input model is not modified
    assert input_model.meta.cal_step.clean_flicker_noise is None

    input_model.close()
    cleaned.close()


def test_save_mask(tmp_path):
    input_model = make_small_ramp_model()
    input_model.meta.filename = "test_jump.fits"

    output_dir = str(tmp_path)
    cleaned = CleanFlickerNoiseStep.call(
        input_model, skip=False, save_mask=True, output_dir=output_dir
    )

    mask_file = str(tmp_path / "test_mask.fits")
    assert os.path.isfile(mask_file)

    mask_model = datamodels.open(mask_file)
    assert isinstance(mask_model, datamodels.ImageModel)

    input_model.close()
    cleaned.close()
    mask_model.close()


def test_save_background(tmp_path):
    input_model = make_small_ramp_model()
    input_model.meta.filename = "test_jump.fits"

    output_dir = str(tmp_path)
    cleaned = CleanFlickerNoiseStep.call(
        input_model, skip=False, save_background=True, output_dir=output_dir
    )

    bkg_file = str(tmp_path / "test_flicker_bkg.fits")
    assert os.path.isfile(bkg_file)

    bkg_model = datamodels.open(bkg_file)
    assert isinstance(bkg_model, datamodels.RampModel)

    input_model.close()
    cleaned.close()
    bkg_model.close()


def test_save_noise(tmp_path):
    input_model = make_small_ramp_model()
    input_model.meta.filename = "test_jump.fits"

    output_dir = str(tmp_path)
    cleaned = CleanFlickerNoiseStep.call(
        input_model, skip=False, save_noise=True, output_dir=output_dir
    )

    noise_file = str(tmp_path / "test_flicker_noise.fits")
    assert os.path.isfile(noise_file)

    noise_model = datamodels.open(noise_file)
    assert isinstance(noise_model, datamodels.RampModel)

    input_model.close()
    cleaned.close()
    noise_model.close()


def test_apply_flat(log_watcher):
    input_model = make_small_ramp_model()

    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise_step", message="Using FLAT")
    cleaned = CleanFlickerNoiseStep.call(input_model, skip=False, apply_flat_field=True)
    watcher.assert_seen()

    # Flat file was used, but flat_field step was not applied
    assert cleaned.meta.ref_file.flat.name is not None
    assert cleaned.meta.cal_step.flat_field is None

    input_model.close()
    cleaned.close()


def test_apply_flat_not_available(log_watcher):
    input_model = make_nirspec_fs_model()

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise_step",
        message="Flat correction is not available",
    )
    cleaned = CleanFlickerNoiseStep.call(input_model, skip=False, apply_flat_field=True)
    watcher.assert_seen()

    # Flat file was not used but step proceeded
    assert cleaned.meta.ref_file.flat.name == "N/A"

    input_model.close()
    cleaned.close()


def test_autoparam_niriss_image(caplog):
    input_model = make_niriss_rate_model()
    cleaned = CleanFlickerNoiseStep.call(input_model, autoparam=True)
    assert "Auto parameters set for NIS_IMAGE" in caplog.text
    assert "apply_flat_field: True" in caplog.text
    assert "Using FLAT reference file" in caplog.text
    assert cleaned.meta.ref_file.flat.name is not None
    assert cleaned.meta.ref_file.flat.name != "N/A"
    assert "background_method: median" in caplog.text


def test_autoparam_nircam_image(caplog):
    input_model = make_nircam_rate_model()
    cleaned = CleanFlickerNoiseStep.call(input_model, autoparam=True)
    assert "Auto parameters set for NRC_IMAGE" in caplog.text
    assert "apply_flat_field: True" in caplog.text
    assert "Using FLAT reference file" in caplog.text
    assert cleaned.meta.ref_file.flat.name is not None
    assert cleaned.meta.ref_file.flat.name != "N/A"
    assert "background_method: median" in caplog.text


def test_autoparam_not_available(caplog):
    input_model = make_small_ramp_model()
    CleanFlickerNoiseStep.call(input_model, autoparam=True)
    assert "Auto parameters are not available for exposure type MIR_IMAGE" in caplog.text
    assert "Using input parameters as provided" in caplog.text


def test_autoparam_failed(caplog, monkeypatch):
    # Mock a failure in autoparam: return None
    monkeypatch.setattr(autoparam, "niriss_image_parameters", lambda *args: None)

    input_model = make_niriss_rate_model()
    CleanFlickerNoiseStep.call(input_model, autoparam=True, override_flat="N/A")
    assert "Auto parameter setting failed" in caplog.text
    assert "Using input parameters as provided" in caplog.text
    assert "apply_flat_field: True" not in caplog.text


def test_output_is_not_input():
    input_model = make_small_ramp_model()
    cleaned = CleanFlickerNoiseStep.call(input_model)

    # successful completion
    assert cleaned.meta.cal_step.clean_flicker_noise == "COMPLETE"

    # make sure input is not modified
    assert cleaned is not input_model
    assert input_model.meta.cal_step.clean_flicker_noise is None

    input_model.close()
    cleaned.close()
