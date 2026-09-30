import numpy as np
from crds.core.exceptions import CrdsLookupError
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.pfpc.pfpc_step import PFPCStep
from jwst.stpipe import Pipeline


def test_step_one_file(mrs_dith1_ch12_medium, mrs_pfpc_model):
    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc=mrs_pfpc_model)

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


def test_step_two_bands(mrs_dith1_ch12_medium, mrs_dith1_ch34_short, mrs_pfpc_model):
    input_models = [mrs_dith1_ch12_medium, mrs_dith1_ch34_short]
    result = PFPCStep.call(input_models, override_pfpc=mrs_pfpc_model)

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


def test_step_two_dithers(mrs_dith1_ch12_medium, mrs_dith2_ch12_medium, mrs_pfpc_model):
    input_models = [mrs_dith1_ch12_medium, mrs_dith2_ch12_medium]
    result = PFPCStep.call(input_models, override_pfpc=mrs_pfpc_model)

    # input is unchanged
    assert mrs_dith1_ch12_medium.meta.cal_step.pfpc is None
    assert mrs_dith2_ch12_medium.meta.cal_step.pfpc is None

    # output is a container with 2 spectra, with dithers averaged by band
    assert isinstance(result, ModelContainer)
    assert len(result) == 2
    for i, spec in enumerate(result):
        assert isinstance(spec, datamodels.MRSMultiSpecModel)
        assert spec.meta.cal_step.pfpc == "COMPLETE"
        assert len(spec.spec) == 1

        # Spectra were extracted from channel 1 and 2 respectively
        assert spec.meta.instrument.channel == str(i + 1)


def test_step_no_pfpc_file(caplog, mrs_dith1_ch12_medium):
    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc="N/A")

    assert "No PFPC reference file found" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_pfpc_lookup_error(caplog, monkeypatch, mrs_dith1_ch12_medium):
    step = PFPCStep()

    def mock_lookup(*args):
        raise CrdsLookupError("No PFPC")

    monkeypatch.setattr(step, "get_reference_file", mock_lookup)

    result = step.run(mrs_dith1_ch12_medium)
    assert "No PFPC reference file found" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_no_matching_correction(caplog, mrs_pfpc_model):
    input_model = datamodels.ImageModel()
    result = PFPCStep.call(input_model, override_pfpc=mrs_pfpc_model)

    assert "No matching PFPC correction" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_not_point_source(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    input_model = mrs_dith1_ch12_medium.copy()
    input_model.meta.target.source_type = "EXTENDED"
    result = PFPCStep.call(input_model, override_pfpc=mrs_pfpc_model)

    assert "Not a point source" in caplog.text
    assert "No correctable spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_no_ta(caplog, monkeypatch, mrs_dith1_ch12_medium, mrs_pfpc_model):
    # TA check not yet implemented: monkeypatch it to fail
    from jwst.pfpc import pfpc as pf

    monkeypatch.setattr(pf, "_ta_performed", lambda x: False)

    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc=mrs_pfpc_model)

    assert "TA was not performed" in caplog.text
    assert "No correctable spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_correction_invalid(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    # modify the PFPC file to trigger failure in correction application

    # no matching correction for channel 1
    idx = mrs_pfpc_model.pfpc_table["channel"] == "1"
    mrs_pfpc_model.pfpc_table["channel"][idx] = "5"

    # all corrections are NaN
    mrs_pfpc_model.pfpc_table["correction"] *= np.nan

    result = PFPCStep.call(mrs_dith1_ch12_medium, override_pfpc=mrs_pfpc_model)

    assert "No matching correction found for test_mrs_ch1" in caplog.text
    assert "No valid correction for test_mrs_ch2" in caplog.text
    assert "No corrected spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_step_combine_invalid(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    # modify the input model so spectra can't be combined
    input_model = mrs_dith1_ch12_medium.copy()
    input_model.meta.exposure.exposure_time = None

    result = PFPCStep.call(input_model, override_pfpc=mrs_pfpc_model)

    assert "No valid spectra were created" in caplog.text

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0


def test_fail_in_pipeline_context(caplog, mrs_dith1_ch12_medium):
    input_copy = mrs_dith1_ch12_medium.copy()
    test_pipe = Pipeline()
    test_step = PFPCStep()
    test_step.parent = test_pipe
    test_step.skip = False

    # override the pfpc file so the step will fail
    test_step.override_pfpc = "N/A"
    result = test_step.run(input_copy)

    # returned result is an empty container
    assert isinstance(result, ModelContainer)
    assert len(result) == 0

    # input has failed status
    assert input_copy.meta.cal_step.pfpc == "FAILED"


def test_output_file(tmp_path, caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    result = PFPCStep.call(
        mrs_dith1_ch12_medium,
        override_pfpc=mrs_pfpc_model,
        save_results=True,
        output_file="test_output_file",
        output_dir=str(tmp_path),
    )
    assert len(result) == 2
    for i, spec in enumerate(result):
        expected = f"test_output_file_ch{i + 1}-medium_pfpc.fits"
        assert spec.meta.filename == expected
        assert (tmp_path / expected).exists()


def test_nirspec_ifu_one_file(caplog, nrs_ifu_dith1, nrs_ifu_pfpc_model):
    result = PFPCStep.call(nrs_ifu_dith1, override_pfpc=nrs_ifu_pfpc_model)

    # input is unchanged
    assert nrs_ifu_dith1.meta.cal_step.pfpc is None

    # output is a container with 1 spectrum
    assert isinstance(result, ModelContainer)
    assert len(result) == 1

    spec = result[0]
    assert isinstance(spec, datamodels.MultiSpecModel)
    assert len(spec.spec) == 1
    assert spec.meta.cal_step.pfpc == "COMPLETE"

    # Output is named with the filter and grating
    assert spec.meta.filename == "test_nrs_ifu_prism-clear"

    # cube_build, extract_1d completed, spectral_leak not attempted
    assert spec.meta.cal_step.cube_build == "COMPLETE"
    assert spec.meta.cal_step.extract_1d == "COMPLETE"
    assert spec.meta.cal_step.spectral_leak is None
