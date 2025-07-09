import json
import logging
import os
from os.path import abspath

import asdf
import pytest
from astropy.extern.configobj.configobj import ConfigObj
from astropy.utils.data import get_pkg_data_filename
from crds.core.exceptions import CrdsLookupError
from stpipe import crds_client
from stpipe import cmdline
from stpipe.config import StepConfig
from stpipe.config_parser import ValidationError
from stdatamodels.jwst import datamodels

from jwst import __version__ as jwst_version
from jwst.white_light import WhiteLightStep
from jwst.tests.helpers import LogWatcher
from jwst.datamodels import ModelLibrary, ModelContainer

from jwst.stpipe import Step
from jwst.stpipe.tests.steps import (
    EmptyPipeline,
    MakeListPipeline,
    MakeListStep,
    ProperPipeline,
    AnotherDummyStep,
)
from jwst.stpipe.tests.steps import OptionalRefTypeStep, SavePipeline

WHITELIGHTSTEP_CRDS_MIRI_PARS = {
    "max_wavelength": 12.0,
    "min_wavelength": 5.0,
    "class": "jwst.white_light.white_light_step.WhiteLightStep",
    "name": "whitelight",
}


CRDS_ERROR_STRING = "PARS-WITHDEFAULTSSTEP: No parameters found"


@pytest.mark.filterwarnings("ignore::ResourceWarning")
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


@pytest.mark.filterwarnings("ignore::ResourceWarning")
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

    data = single_member_asn
    if on_disk_status is not None:
        with ModelLibrary(single_member_asn, on_disk=on_disk_status) as model_library:
            data = model_library

    pars = WhiteLightStep.get_config_from_reference(data)
    assert pars == WHITELIGHTSTEP_CRDS_MIRI_PARS


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


@pytest.mark.parametrize(
    "cfg_file, expected_reftype",
    [
        ("local_class.cfg", "pars-dummystep"),
        ("jwst_generic_pars-makeliststep_0002.asdf", "pars-makeliststep"),
    ],
)
def test_reftype(cfg_file, expected_reftype):
    """Test that reftype is produced as expected"""
    step = Step.from_config_file(
        get_pkg_data_filename(f"steps/{cfg_file}", package="jwst.stpipe.tests")
    )
    assert step.__class__.get_config_reftype() == expected_reftype
    assert step.get_config_reftype() == expected_reftype


def test_saving_pars(tmp_path):
    """Save the step parameters from the commandline"""
    cfg_path = get_pkg_data_filename(
        "steps/jwst_generic_pars-makeliststep_0002.asdf", package="jwst.stpipe.tests"
    )
    saved_path = str(tmp_path / "savepars.asdf")
    Step.from_cmdline([cfg_path, "--save-parameters", str(saved_path)])
    assert os.path.exists(saved_path)

    with asdf.open(
        get_pkg_data_filename(
            "steps/jwst_generic_pars-makeliststep_0002.asdf", package="jwst.stpipe.tests"
        )
    ) as af:
        original_config = StepConfig.from_asdf(af)
        original_config.parameters["par3"] = False

    with asdf.open(str(saved_path)) as af:
        config = StepConfig.from_asdf(af)
        assert config.parameters == original_config.parameters


@pytest.mark.parametrize(
    "step_obj, expected",
    [
        (
            MakeListStep(par1=0.0, par2="from args"),
            StepConfig(
                "jwst.stpipe.tests.steps.MakeListStep",
                "MakeListStep",
                {
                    "pre_hooks": [],
                    "post_hooks": [],
                    "output_ext": ".fits",
                    "output_use_model": False,
                    "output_use_index": True,
                    "save_results": False,
                    "skip": False,
                    "search_output_file": True,
                    "input_dir": "",
                    "par1": 0.0,
                    "par2": "from args",
                    "par3": False,
                },
                [],
            ),
        ),
    ],
)
def test_export_config(step_obj, expected, tmp_path):
    """Test retrieving of configuration parameters"""
    config_path = tmp_path / "config.asdf"
    step_obj.export_config(config_path)

    with asdf.open(config_path) as af:
        # StepConfig has an __eq__ implementation but we can't use it
        # due to differences between asdf 2.7 and 2.8 in serializing None
        # values.  This can be simplified once the minimum asdf requirement
        # is changed to >= 2.8.
        # assert StepConfig.from_asdf(af) == expected
        config = StepConfig.from_asdf(af)
        assert config.class_name == expected.class_name
        assert config.name == expected.name
        assert config.steps == expected.steps
        parameters = set(expected.parameters.keys()).union(set(config.parameters.keys()))
        for parameter in parameters:
            assert config.parameters.get(parameter) == expected.parameters.get(parameter)


@pytest.mark.parametrize(
    "step_obj, full_spec, expected",
    [
        # #############################################
        # Test with `full_spec = True`
        # #############################################
        # Mix of required and optional parameters
        (
            MakeListStep(par1=0.0, par2="from args"),
            True,
            {
                "pre_hooks": [],
                "post_hooks": [],
                "output_file": None,
                "output_dir": None,
                "output_ext": ".fits",
                "output_use_model": False,
                "output_use_index": True,
                "save_results": False,
                "skip": False,
                "suffix": None,
                "search_output_file": True,
                "input_dir": "",
                "par1": 0.0,
                "par2": "from args",
                "par3": False,
            },
        ),
        #
        # All parameters set
        #
        (
            MakeListPipeline(
                par1="Instantiated", steps={"make_list": {"par1": 0.0, "par2": "sub-instantiated"}}
            ),
            True,
            {
                "pre_hooks": [],
                "post_hooks": [],
                "output_file": None,
                "output_dir": None,
                "output_ext": ".fits",
                "output_use_model": False,
                "output_use_index": True,
                "save_results": False,
                "skip": False,
                "suffix": None,
                "search_output_file": True,
                "input_dir": "",
                "par1": "Instantiated",
                "steps": {
                    "make_list": {
                        "pre_hooks": [],
                        "post_hooks": [],
                        "output_file": None,
                        "output_dir": None,
                        "output_ext": ".fits",
                        "output_use_model": False,
                        "output_use_index": True,
                        "save_results": False,
                        "skip": False,
                        "suffix": None,
                        "search_output_file": True,
                        "input_dir": "",
                        "par1": 0.0,
                        "par2": "sub-instantiated",
                        "par3": False,
                    }
                },
            },
        ),
        #
        # Pipeline without any sub-steps
        #
        (
            EmptyPipeline(par1="Instantiated"),
            True,
            {
                "pre_hooks": [],
                "post_hooks": [],
                "output_file": None,
                "output_dir": None,
                "output_ext": ".fits",
                "output_use_model": False,
                "output_use_index": True,
                "save_results": False,
                "skip": False,
                "suffix": None,
                "search_output_file": True,
                "input_dir": "",
                "par1": "Instantiated",
                "steps": {},
            },
        ),
        # ######################################
        # Test with `full_spec=False`
        # ######################################
        # Mix of required and optional parameters
        (
            MakeListStep(par1=0.0, par2="from args"),
            False,
            {"par1": 0.0, "par2": "from args", "par3": False},
        ),
        # Pipeline with all parameters set
        (
            MakeListPipeline(
                par1="Instantiated", steps={"make_list": {"par1": 0.0, "par2": "sub-instantiated"}}
            ),
            False,
            {
                "par1": "Instantiated",
                "steps": {"make_list": {"par1": 0.0, "par2": "sub-instantiated", "par3": False}},
            },
        ),
        # Pipeline without any sub-steps
        (EmptyPipeline(par1="Instantiated"), False, {"par1": "Instantiated", "steps": {}}),
    ],
)
def test_getpars(step_obj, full_spec, expected):
    """Test retrieving of configuration parameters"""
    assert step_obj.get_pars(full_spec=full_spec) == expected


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


@pytest.mark.parametrize(
    "command_line_pars, command_line_config_pars, reference_pars, expected_pars",
    [
        # If nothing else is present, we should use the spec defaults
        (
            None,
            None,
            None,
            {
                "par1": "default par1 value",
                "par2": "default par2 value",
                "par3": "default par3 value",
                "par4": "default par4 value",
            },
        ),
        # Reference file pars > spec defaults
        (
            None,
            None,
            {
                "par1": "reference par1 value",
                "par2": "reference par2 value",
                "par3": "reference par3 value",
                "par4": "reference par4 value",
            },
            {
                "par1": "reference par1 value",
                "par2": "reference par2 value",
                "par3": "reference par3 value",
                "par4": "reference par4 value",
            },
        ),
        # Command line config pars > reference pars
        (
            None,
            {
                "par1": "config par1 value",
                "par2": "config par2 value",
                "par3": "config par3 value",
                "par4": "config par4 value",
            },
            {
                "par1": "reference par1 value",
                "par2": "reference par2 value",
                "par3": "reference par3 value",
                "par4": "reference par4 value",
            },
            {
                "par1": "config par1 value",
                "par2": "config par2 value",
                "par3": "config par3 value",
                "par4": "config par4 value",
            },
        ),
        # Command line override pars > all other pars
        (
            {
                "par1": "override par1 value",
                "par2": "override par2 value",
                "par3": "override par3 value",
                "par4": "override par4 value",
            },
            {
                "par1": "config par1 value",
                "par2": "config par2 value",
                "par3": "config par3 value",
                "par4": "config par4 value",
            },
            {
                "par1": "reference par1 value",
                "par2": "reference par2 value",
                "par3": "reference par3 value",
                "par4": "reference par4 value",
            },
            {
                "par1": "override par1 value",
                "par2": "override par2 value",
                "par3": "override par3 value",
                "par4": "override par4 value",
            },
        ),
        # Test complex merging of parameters (one parameter per source)
        (
            {"par1": "override par1 value"},
            {"par2": "config par2 value"},
            {"par3": "reference par3 value"},
            {
                "par1": "override par1 value",
                "par2": "config par2 value",
                "par3": "reference par3 value",
                "par4": "default par4 value",
            },
        ),
        # Test complex merging of parameters (more parameters specified on the lower-precedence sources)
        (
            {"par1": "override par1 value"},
            {"par1": "config par1 value", "par2": "config par2 value"},
            {
                "par1": "reference par1 value",
                "par2": "reference par2 value",
                "par3": "reference par3 value",
            },
            {
                "par1": "override par1 value",
                "par2": "config par2 value",
                "par3": "reference par3 value",
                "par4": "default par4 value",
            },
        ),
    ],
)
def test_step_from_commandline_par_precedence(
    command_line_pars,
    command_line_config_pars,
    reference_pars,
    expected_pars,
    tmp_path,
    monkeypatch,
):
    args = []

    class_name = "jwst.stpipe.tests.steps.WithDefaultsStep"
    config_name = "WithDefaultsStep"
    reference_type = f"pars-{config_name.lower()}"
    input_path = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")

    if command_line_config_pars:
        command_line_config_path = tmp_path / "with_defaults_step.cfg"
        config = ConfigObj(str(command_line_config_path))
        config["class"] = class_name
        config["name"] = config_name
        for key, value in command_line_config_pars.items():
            config[key] = value
        config.write()
        args.append(str(command_line_config_path.absolute()))
    else:
        args.append(class_name)

    args.append(input_path)

    if command_line_pars:
        for key, value in command_line_pars.items():
            args.append(f"--{key}={value}")

    reference_file_map = {}
    if reference_pars:
        reference_path = tmp_path / f"{reference_type}.asdf"
        reference_config = StepConfig(class_name, config_name, reference_pars, [])
        with reference_config.to_asdf() as af:
            af.write_to(reference_path)

        reference_file_map[reference_type] = str(reference_path)

    def mock_get_reference_file(dataset, reference_file_type, observatory=None, asn_exptypes=None):
        if reference_file_type in reference_file_map:
            return reference_file_map[reference_file_type]
        else:
            raise CrdsLookupError(
                f"Error determining best reference for '{reference_file_type}'  = \
  Unknown reference type '{reference_file_type}'"
            )

    monkeypatch.setattr(crds_client, "get_reference_file", mock_get_reference_file)

    step = Step.from_cmdline(args)

    for key, value in expected_pars.items():
        assert getattr(step, key) == value


def test_step_with_local_class():
    step_fn = get_pkg_data_filename("steps/local_class.cfg", package="jwst.stpipe.tests")
    step = Step.from_config_file(step_fn)

    step.run(datamodels.ImageModel((2, 2)))


def test_extra_parameter():
    with pytest.raises(ValidationError):
        AnotherDummyStep("SomeOtherStepOriginal", par5="foo")


def test_crds_override():
    ff_name = get_pkg_data_filename("data/flat.fits", package="jwst.stpipe.tests")
    step = AnotherDummyStep(
        "SomeOtherStepOriginal", par1=42.0, par2="abc def", override_flat_field=ff_name
    )

    fd = step.get_reference_file(datamodels.open(), "flat_field")
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


def test_dunder_call_error():
    pipeline = EmptyPipeline()
    with pytest.raises(TypeError, match="not callable"):
        pipeline(None)
