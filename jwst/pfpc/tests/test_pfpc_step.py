import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.pfpc.pfpc_step import PFPCStep
from jwst.pfpc.tests import helpers


@pytest.fixture(scope="module")
def mrs_dith1_ch12_medium():
    model = helpers.dithered_miri_mrs_model(detector="MIRIFUSHORT", channel="12", band="MEDIUM")
    yield model
    model.close()


@pytest.fixture(scope="module")
def mrs_dith1_ch34_short():
    model = helpers.dithered_miri_mrs_model(detector="MIRIFULONG", channel="34", band="SHORT")
    yield model
    model.close()


@pytest.fixture(scope="module")
def pfpc_model():
    pfpc = helpers.miri_mrs_pfpc_model()
    yield pfpc
    pfpc.close()


def test_step_one_file(mrs_dith1_ch12_medium, pfpc_model):
    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc=pfpc_model)

    # input is unchanged
    assert mrs_dith1_ch12_medium.meta.cal_step.pfpc is None

    # output is a container with 2 spectra
    assert isinstance(result, ModelContainer)
    assert len(result) == 2
    for i, spec in enumerate(result):
        assert isinstance(spec, datamodels.MRSMultiSpecModel)
        assert len(spec.spec) == 1
        assert spec.meta.cal_step.pfpc == "COMPLETE"

        # Spectra were extracted from channel 1 and 2 respectively
        assert spec.meta.instrument.channel == str(i + 1)
        assert spec.meta.filename == f"test_mrs_ch{i + 1}-medium"

        # cube_build, extract_1d completed, spectral_leak skipped
        assert spec.meta.cal_step.cube_build == "COMPLETE"
        assert spec.meta.cal_step.extract_1d == "COMPLETE"
        assert spec.meta.cal_step.spectral_leak == "SKIPPED"


def test_step_two_files(mrs_dith1_ch12_medium, mrs_dith1_ch34_short, pfpc_model):
    input_models = [mrs_dith1_ch12_medium, mrs_dith1_ch34_short]
    result = PFPCStep.call(input_models, override_pfpc=pfpc_model)

    # input is unchanged
    assert mrs_dith1_ch12_medium.meta.cal_step.pfpc is None
    assert mrs_dith1_ch34_short.meta.cal_step.pfpc is None

    # output is a container with 4 spectra
    assert isinstance(result, ModelContainer)
    assert len(result) == 4
    for i, spec in enumerate(result):
        assert isinstance(spec, datamodels.MRSMultiSpecModel)
        assert len(spec.spec) == 1
        assert spec.meta.cal_step.pfpc == "COMPLETE"

        # Spectra were extracted from channels 1-4 respectively
        assert spec.meta.instrument.channel == str(i + 1)
        if i < 2:
            assert spec.meta.filename == f"test_mrs_ch{i + 1}-medium"
        else:
            assert spec.meta.filename == f"test_mrs_ch{i + 1}-short"

        # cube_build, extract_1d completed
        assert spec.meta.cal_step.cube_build == "COMPLETE"
        assert spec.meta.cal_step.extract_1d == "COMPLETE"

        # spectral_leak also completed for channel 3 since the input contains ch1b and ch3a
        if spec.meta.instrument.channel == "3":
            assert spec.meta.cal_step.spectral_leak == "COMPLETE"
        else:
            # For the other channels, the status is not updated
            assert spec.meta.cal_step.spectral_leak is None


def test_step_no_pfpc_file(caplog, mrs_dith1_ch12_medium):
    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc="N/A")

    assert "No PFPC reference file found" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_no_matching_correction(caplog, pfpc_model):
    input_model = datamodels.ImageModel()
    result = PFPCStep.call(input_model, override_pfpc=pfpc_model)

    assert "No matching PFPC correction" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_not_point_source(caplog, mrs_dith1_ch12_medium, pfpc_model):
    input_model = mrs_dith1_ch12_medium.copy()
    input_model.meta.target.source_type = "EXTENDED"
    result = PFPCStep.call(input_model, override_pfpc=pfpc_model)

    assert "Not a point source" in caplog.text
    assert "No correctable spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_correction_invalid(caplog, mrs_dith1_ch12_medium, pfpc_model):
    # modify the PFPC file to trigger failure in correction application
    pfpc_copy = pfpc_model.copy()

    # no matching correction for channel 1
    idx = pfpc_copy.pfpc_table["channel"] == "1"
    pfpc_copy.pfpc_table["channel"][idx] = "5"

    # all corrections are NaN
    pfpc_copy.pfpc_table["correction"] *= np.nan

    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc=pfpc_copy)

    assert "No matching correction found for test_mrs_ch1" in caplog.text
    assert "No valid correction for test_mrs_ch2" in caplog.text
    assert "No corrected spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_combine_invalid(caplog, mrs_dith1_ch12_medium, pfpc_model):
    # modify the input model so spectra can't be combined
    input_model = mrs_dith1_ch12_medium.copy()
    input_model.meta.exposure.exposure_time = None

    result = PFPCStep.call(input_model, override_pfpc=pfpc_model)

    assert "No valid spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0
