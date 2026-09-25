import os

import pytest
from stdatamodels.jwst import datamodels

from jwst.clean_flicker_noise.tests.helpers import (
    make_nirspec_fs_model,
    make_nrs_fs_full_ramp,
    make_small_ramp_model,
)
from jwst.picture_frame import PictureFrameStep
from jwst.picture_frame.tests.helpers import picture_frame_model
from jwst.pipeline import Detector1Pipeline


@pytest.fixture()
def picture_frame_file(tmp_path):
    model = picture_frame_model()
    pctfrm_file = str(tmp_path / "picture_frame.fits")
    model.save(pctfrm_file)
    model.close()
    return pctfrm_file


@pytest.mark.parametrize("skip", [True, False])
def test_output_type(skip, picture_frame_file):
    input_model = make_nrs_fs_full_ramp()

    # Default is that step is skipped
    step_instance = PictureFrameStep()
    assert step_instance.skip == True

    # Run the step with skip option
    step_instance.skip = skip
    step_instance.override_pictureframe = picture_frame_file
    cleaned = step_instance.run(input_model)

    # Since it's run standalone, the step is completed either way
    assert isinstance(cleaned, datamodels.RampModel)
    assert cleaned.meta.cal_step.picture_frame == "COMPLETE"

    # Input model is not modified
    assert cleaned is not input_model
    assert input_model.meta.cal_step.picture_frame is None

    input_model.close()
    cleaned.close()


@pytest.mark.parametrize("skip", [True, False])
def test_run_in_pipeline(skip, picture_frame_file):
    input_model = make_nrs_fs_full_ramp()

    # Default is to skip the step
    pipeline_instance = Detector1Pipeline()
    assert pipeline_instance.steps["picture_frame"]["skip"] == True

    # Run the pipeline, skipping most steps for processing time
    pipeline_instance.saturation.skip = True
    pipeline_instance.ipc.skip = True
    pipeline_instance.superbias.skip = True
    pipeline_instance.refpix.skip = True
    pipeline_instance.linearity.skip = True
    pipeline_instance.persistence.skip = True
    pipeline_instance.dark_current.skip = True
    pipeline_instance.jump.skip = True
    pipeline_instance.clean_flicker_noise.skip = True
    pipeline_instance.ramp_fit.skip = True
    pipeline_instance.gain_scale.skip = True

    # Run the picture frame step if needed
    pipeline_instance.picture_frame.skip = skip
    pipeline_instance.picture_frame.override_pictureframe = picture_frame_file

    cleaned = pipeline_instance.run(input_model)

    if skip:
        assert cleaned.meta.cal_step.picture_frame == "SKIPPED"
    else:
        assert cleaned.meta.cal_step.picture_frame == "COMPLETE"

    # Either way, input model should not be modified
    assert cleaned is not input_model
    assert input_model.meta.cal_step.picture_frame is None

    input_model.close()
    cleaned.close()


def test_save_mask(tmp_path, picture_frame_file):
    input_model = make_nirspec_fs_model()
    input_model.meta.filename = "test_rate.fits"

    output_dir = str(tmp_path)
    cleaned = PictureFrameStep.call(
        input_model,
        skip=False,
        save_mask=True,
        output_dir=output_dir,
        override_pictureframe=picture_frame_file,
    )

    mask_file = str(tmp_path / "test_pctfrm_mask.fits")
    assert os.path.isfile(mask_file)

    mask_model = datamodels.open(mask_file)
    assert isinstance(mask_model, datamodels.ImageModel)

    input_model.close()
    cleaned.close()
    mask_model.close()


def test_save_correction(tmp_path, picture_frame_file):
    input_model = make_nirspec_fs_model()
    input_model.meta.filename = "test_rate.fits"

    output_dir = str(tmp_path)
    cleaned = PictureFrameStep.call(
        input_model,
        skip=False,
        save_correction=True,
        output_dir=output_dir,
        override_pictureframe=picture_frame_file,
    )

    corr_file = str(tmp_path / "test_pctfrm_correction.fits")
    assert os.path.isfile(corr_file)

    corr_model = datamodels.open(corr_file)
    assert isinstance(corr_model, datamodels.ImageModel)

    input_model.close()
    cleaned.close()
    corr_model.close()


def test_unsupported_mode(caplog):
    # Make a MIRI model - it will never have a picture frame file in CRDS
    input_model = make_small_ramp_model()
    cleaned = PictureFrameStep.call(input_model)

    assert "Picture frame correction is not available" in caplog.text
    assert cleaned.meta.cal_step.picture_frame == "SKIPPED"

    # Input is not modified
    assert cleaned is not input_model
    assert input_model.meta.cal_step.picture_frame is None
