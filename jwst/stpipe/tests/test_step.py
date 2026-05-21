import logging
from os.path import abspath

import pytest
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels
from stpipe.config_parser import ValidationError

from jwst import __version__ as jwst_version
from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe import Pipeline, Step
from jwst.stpipe.tests.steps import (
    AnotherDummyStep,
    EmptyPipeline,
    OptionalRefTypeStep,
    ProperPipeline,
    SavePipeline,
)
from jwst.tests.helpers import LogWatcher
from jwst.white_light import WhiteLightStep

WHITELIGHTSTEP_CRDS_MIRI_PARS = {
    "max_wavelength": 12.0,
    "min_wavelength": 5.0,
    "class": "jwst.white_light.white_light_step.WhiteLightStep",
    "name": "whitelight",
}


CRDS_ERROR_STRING = "PARS-WITHDEFAULTSSTEP: No parameters found"


@pytest.mark.parametrize(
    "arg, env_set, expected_fn",
    [
        ("--disable-crds-steppars", None, lambda stream: CRDS_ERROR_STRING not in stream),
        ("--verbose", None, lambda stream: CRDS_ERROR_STRING in stream),
        ("--verbose", "true", lambda stream: CRDS_ERROR_STRING not in stream),
        ("--verbose", "True", lambda stream: CRDS_ERROR_STRING not in stream),
        ("--verbose", "t", lambda stream: CRDS_ERROR_STRING not in stream),
    ],
)
def test_disable_crds_steppars_cmdline(capsys, arg, env_set, expected_fn, monkeypatch):
    """Test setting of disable_crds_steppars"""
    data_path = get_pkg_data_filename("data/miri_data.fits", package="jwst.stpipe.tests")

    if env_set:
        monkeypatch.setenv("STPIPE_DISABLE_CRDS_STEPPARS", env_set)

    Step.from_cmdline(["jwst.stpipe.tests.steps.WithDefaultsStep", data_path, arg])

    captured = capsys.readouterr()
    assert expected_fn(captured.err)


def test_parameters_from_crds_open_model():
    """Test retrieval of parameters from CRDS"""
    with datamodels.open(
        get_pkg_data_filename("data/miri_data.fits", package="jwst.stpipe.tests")
    ) as data:
        pars = WhiteLightStep.get_config_from_reference(data)
    assert pars == WHITELIGHTSTEP_CRDS_MIRI_PARS


def test_parameters_from_crds_filename(monkeypatch):
    """
    Test retrieval of parameters from CRDS from single-model filename.

    Ensures datamodels.open() is not called.
    Similar tests of read_metadata() in stdatamodels ensure that the input file's `data`
    attribute is never read either.
    """

    def throw_error(self):
        raise Exception()  # noqa: TRY002

    monkeypatch.setattr(datamodels, "open", throw_error)

    fname = get_pkg_data_filename("data/miri_data.fits", package="jwst.stpipe.tests")
    pars = WhiteLightStep.get_config_from_reference(fname)
    assert pars == WHITELIGHTSTEP_CRDS_MIRI_PARS


@pytest.mark.parametrize("on_disk_status", [None, True, False])
def test_parameters_from_crds_association(on_disk_status, monkeypatch):
    """
    Test retrieval of parameters from CRDS from an association or library.

    If on_disk_status is not None, the association is first opened as a library
    then passed to get_config_from_reference with on_disk set to on_disk_status.
    Otherwise, the asn file is passed directly to get_config_from_reference.

    The datamodel should never be opened in any of the three cases, even if
    on_disk=False, because the model has not yet been borrowed from the library.
    """

    def throw_error(self):
        raise Exception()  # noqa: TRY002

    monkeypatch.setattr(datamodels, "open", throw_error)

    single_member_asn = get_pkg_data_filename(
        "data/single_member_miri_asn.json", package="jwst.stpipe.tests"
    )

    if on_disk_status is not None:
        with ModelLibrary(single_member_asn, on_disk=on_disk_status) as model_library:
            pars = WhiteLightStep.get_config_from_reference(model_library)
    else:
        pars = WhiteLightStep.get_config_from_reference(single_member_asn)

    assert pars == WHITELIGHTSTEP_CRDS_MIRI_PARS

    if on_disk_status is True:
        model_library._temp_dir.cleanup()  # avoid ResourceWarning


@pytest.mark.parametrize("is_list", [True, False])
def test_parameters_from_crds_listlike(is_list):
    """Test retrieval of parameters from CRDS from a list of open models or ModelContainer"""
    single_member_asn = get_pkg_data_filename(
        "data/single_member_miri_asn.json", package="jwst.stpipe.tests"
    )
    data = ModelContainer(single_member_asn)
    if is_list:
        data = data._models
    pars = WhiteLightStep.get_config_from_reference(data)
    assert pars == WHITELIGHTSTEP_CRDS_MIRI_PARS

    for model in data:
        model.close()


@pytest.mark.parametrize("is_list", [True, False])
def test_parameters_from_crds_empty_listlike(is_list):
    """Test failure when attempting retrieval of CRDS parameters from an empty list or container"""
    if is_list:
        data = []
    else:
        data = ModelContainer()
    pars = WhiteLightStep.get_config_from_reference(data)
    assert not len(pars)


def test_parameters_from_crds_bad_type():
    """Test failure when attempting retrieval of CRDS parameters from an invalid type"""
    pars = WhiteLightStep.get_config_from_reference(42)
    assert not len(pars)


def test_parameters_from_crds_bad_meta():
    """Test failure when attempting retrieval of parameters from CRDS with invalid metadata"""
    with datamodels.open(
        get_pkg_data_filename("data/miri_data.fits", package="jwst.stpipe.tests")
    ) as data:
        data.meta.instrument.name = "NIRSPEC"
        pars = WhiteLightStep.get_config_from_reference(data)
    assert not len(pars)


def test_hook():
    """Test the running of hooks"""
    step_fn = get_pkg_data_filename("steps/stepwithmodel_hook.cfg", package="jwst.stpipe.tests")
    step = Step.from_config_file(step_fn)

    model = datamodels.ImageModel()
    result = step.run(model)

    assert result.pre_hook_run
    assert step.pre_hook_run
    assert result.post_hook_run
    assert step.post_hook_run


def test_hook_with_return():
    """Test the running of hooks"""
    step_fn = get_pkg_data_filename(
        "steps/stepwithmodel_hookreturn.cfg", package="jwst.stpipe.tests"
    )
    step = Step.from_config_file(step_fn)

    model = datamodels.ImageModel()
    result = step.run(model)

    assert result == "PostHookWithReturnStep executed"
    assert step.pre_hook_run
    assert step.post_hook_run


def test_step():
    step_fn = get_pkg_data_filename("steps/some_other_step.cfg", package="jwst.stpipe.tests")
    step = Step.from_config_file(step_fn)

    assert isinstance(step, AnotherDummyStep)
    assert step.name == "SomeOtherStepOriginal"
    assert step.par2 == "abc def"

    step.run(1, 2)


def test_step_from_python():
    step = AnotherDummyStep("SomeOtherStepOriginal", par1=42.0, par2="abc def")

    assert step.par1 == 42.0
    assert step.par2 == "abc def"
    assert step.par3 is False

    result = step.run(1, 2)

    assert result == 3


def test_step_from_python_simple():
    result = AnotherDummyStep.call(1, 2, par1=42.0, par2="abc def")

    assert result == 3


def test_step_from_python_simple2():
    step_fn = get_pkg_data_filename("steps/some_other_step.cfg", package="jwst.stpipe.tests")

    result = AnotherDummyStep.call(1, 2, config_file=step_fn)

    assert result == 3


def test_step_from_commandline():
    args = [
        abspath(get_pkg_data_filename("steps/some_other_step.cfg", package="jwst.stpipe.tests")),
        "--par1=58",
        "--par2=hij klm",
    ]

    step = Step.from_cmdline(args)

    assert step.par1 == 58.0
    assert step.par2 == "hij klm"
    assert step.par3 is True

    step.run(1, 2)


def test_step_from_commandline_class():
    args = ["jwst.stpipe.tests.steps.AnotherDummyStep", "--par1=58", "--par2=hij klm"]

    step = Step.from_cmdline(args)

    assert step.par1 == 58.0
    assert step.par2 == "hij klm"
    assert step.par3 is False

    step.run(1, 2)


def test_step_from_commandline_class_alias(mock_stpipe_entry_points):
    args = ["stpipe_dummy", "--par1=58", "--par2=hij klm"]

    step = Step.from_cmdline(args)

    assert isinstance(step, AnotherDummyStep)

    assert step.par1 == 58.0
    assert step.par2 == "hij klm"
    assert step.par3 is False

    step.run(1, 2)


def test_step_from_commandline_config_class_alias(mock_stpipe_entry_points):
    args = [
        abspath(get_pkg_data_filename("steps/dummy_with_alias.cfg", package="jwst.stpipe.tests")),
        "--par1=58",
        "--par2=hij klm",
    ]

    step = Step.from_cmdline(args)

    assert step.par1 == 58.0
    assert step.par2 == "hij klm"
    assert step.par3 is True

    step.run(1, 2)


def test_step_from_commandline_invalid():
    args = ["__foo__"]

    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_from_commandline_invalid2():
    args = ["__foo__.__bar__"]
    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_from_commandline_invalid3():
    args = ["sys.foo"]
    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_from_commandline_invalid4():
    args = ["sys.argv"]
    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_extra_parameter():
    with pytest.raises(ValidationError):
        AnotherDummyStep("SomeOtherStepOriginal", par5="foo")


def test_crds_override():
    ff_name = get_pkg_data_filename("data/flat.fits", package="jwst.stpipe.tests")
    step = AnotherDummyStep(
        "SomeOtherStepOriginal", par1=42.0, par2="abc def", override_flat_field=ff_name
    )

    fd = step.get_reference_file(datamodels.JwstDataModel(), "flat_field")
    assert fd == ff_name


def test_omit_ref_file():
    step = OptionalRefTypeStep(override_to_be_ignored_ref_type="")
    step.process()


def test_search_attr(tmp_path):
    value = str(tmp_path)
    pipeline = SavePipeline("afile.fits", output_dir=value)  # codespell:ignore afile

    assert pipeline.search_attr("output_dir") == value
    assert pipeline.stepwithmodel.search_attr("output_dir") == value
    assert pipeline.search_attr("junk") is None
    assert pipeline.stepwithmodel.search_attr("junk") is None


def test_print_configspec():
    step = Step()
    step.print_configspec()


def test_call_with_config(caplog, tmp_cwd):
    """Test call using a config file with substeps

    In particular, from JP-1482, there was a case where a substep parameter
    was not being overridden. Test for that case.
    """
    cfg = get_pkg_data_filename("data/proper_pipeline.asdf", package="jwst.stpipe.tests")
    model = get_pkg_data_filename("data/flat.fits", package="jwst.stpipe.tests")

    ProperPipeline.call(model, config_file=cfg)

    assert "newpar1" in caplog.text


def test_finalize_logging(monkeypatch):
    """
    Check that the jwst version and crds context are logged
    when a step/pipeline is run.
    """
    pipeline = EmptyPipeline()
    model = datamodels.ImageModel()
    watcher = LogWatcher(f"Results used jwst version: {jwst_version}")
    monkeypatch.setattr(logging.getLogger("jwst.stpipe.core"), "info", watcher)
    pipeline.run(model)
    assert watcher.seen


def test_add_asn_id_to_output_name_from_step():
    asn_id_step = "1234"

    # Bare step and model, no ASN ID set yet
    step = Step()
    model = datamodels.ImageModel()

    # No ASN ID in the default output path
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == "test_cal.fits"

    # If input has no ASN ID, the path is unmodified
    found_asn_id = step.add_asn_id_to_output_name(model)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == "test_cal.fits"
    assert found_asn_id is None

    # Add the ASN ID to the step attributes: it appears in the output name
    step.asn_id = asn_id_step
    found_asn_id = step.add_asn_id_to_output_name(model)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == f"test_{asn_id_step}_cal.fits"
    assert found_asn_id == asn_id_step

    # Reset the ASN ID and update again: it no longer appears in the output name
    step.asn_id = None
    found_asn_id = step.add_asn_id_to_output_name(model)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == "test_cal.fits"
    assert found_asn_id is None


def test_add_asn_id_to_output_name_from_step_parent():
    asn_id_step = "1234"
    asn_id_parent = "2345"

    # Bare step and model, parent has ASN ID set
    pipe = Pipeline()
    pipe.asn_id = asn_id_parent
    step = Step(parent=pipe)
    model = datamodels.ImageModel()

    # Parent ASN ID appears in output name
    found_asn_id = step.add_asn_id_to_output_name(model)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == f"test_{asn_id_parent}_cal.fits"
    assert found_asn_id == asn_id_parent

    # If the step has an asn_id defined, that overrides the parent
    step.asn_id = asn_id_step
    found_asn_id = step.add_asn_id_to_output_name(model)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == f"test_{asn_id_step}_cal.fits"
    assert found_asn_id == asn_id_step


@pytest.mark.parametrize("models", ["container", "library"])
def test_add_asn_id_to_output_name_from_model(models):
    asn_id_step = "1234"
    asn_id_model = "2345"

    # Make input models with ASN ID assigned
    if models == "container":
        input_models = ModelContainer()
        input_models.asn_table["asn_id"] = asn_id_model
    elif models == "library":
        model = datamodels.ImageModel()
        input_models = ModelLibrary([model])
        input_models._asn["asn_id"] = asn_id_model

    # Step with no asn_id attribute set
    step = Step()

    # Add the ASN ID from the model
    found_asn_id = step.add_asn_id_to_output_name(input_models)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == f"test_{asn_id_model}_cal.fits"
    assert found_asn_id == asn_id_model

    # If the step has an asn_id, it still uses the one from the model
    step.asn_id = asn_id_step
    found_asn_id = step.add_asn_id_to_output_name(input_models)
    output_name = step.make_output_path(basepath="test", suffix="cal")
    assert output_name == f"test_{asn_id_model}_cal.fits"
    assert found_asn_id == asn_id_model
