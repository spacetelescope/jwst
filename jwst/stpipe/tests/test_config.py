from datetime import datetime, timedelta, timezone

import pytest
import asdf

from stpipe.config import export_config, StepConfig, _validate_asdf

from jwst.stpipe.tests.steps import MakeListPipeline, WithDefaultsStep


@pytest.fixture
def config():
    s1 = StepConfig("some.StepClass", "Step1", {"woo": "hoo"}, [])
    s2 = StepConfig("some.other.StepClass", "Step2", {"yee": "haw"}, [])

    return StepConfig("some.PipelineClass", "Pipeline", {"foo": "bar", "baz": 42}, [s1, s2])


def test_step_config_equality(config):
    other = StepConfig(config.class_name, config.name, config.parameters, config.steps)
    assert config == other

    other = StepConfig("some.other.pipeline.Name", config.name, config.parameters, config.steps)
    assert config != other

    other = StepConfig(config.class_name, "Hildegaard", config.parameters, config.steps)
    assert config != other

    other = StepConfig(config.class_name, config.name, {"foo": "foz"}, config.steps)
    assert config != other

    other = StepConfig(config.class_name, config.name, config.parameters, [config.steps[0]])
    assert config != other

    assert config != object()


def test_step_config_to_asdf(config):
    # Without metadata
    asdf_file = config.to_asdf()
    assert isinstance(asdf_file, asdf.AsdfFile)
    assert asdf_file["class"] == config.class_name
    assert asdf_file["name"] == config.name
    assert asdf_file["parameters"] == config.parameters
    assert asdf_file["steps"] == [s.to_asdf().tree for s in config.steps]
    assert "meta" not in asdf_file
    # No validation errors:
    with asdf.AsdfFile(asdf_file.tree) as af:
        _validate_asdf(af, "http://stsci.edu/schemas/stpipe/step_config-1.0.0")

    # With metadata
    asdf_file = config.to_asdf(include_metadata=True)
    assert isinstance(asdf_file, asdf.AsdfFile)
    assert asdf_file["class"] == config.class_name
    assert asdf_file["name"] == config.name
    assert asdf_file["parameters"] == config.parameters
    assert asdf_file["steps"] == [s.to_asdf().tree for s in config.steps]
    assert asdf_file["meta"]["author"] == "<SPECIFY>"
    current_time = datetime.now(timezone.utc).replace(tzinfo=None)
    assert (current_time - datetime.fromisoformat(asdf_file["meta"]["date"])) < timedelta(seconds=10)
    assert asdf_file["meta"]["description"] == "Parameters for calibration step some.PipelineClass"
    assert asdf_file["meta"]["instrument"]["name"] == "<SPECIFY>"
    assert asdf_file["meta"]["origin"] == "<SPECIFY>"
    assert asdf_file["meta"]["pedigree"] == "<SPECIFY>"
    assert asdf_file["meta"]["reftype"] == "<SPECIFY>"
    assert asdf_file["meta"]["telescope"] == "<SPECIFY>"
    assert asdf_file["meta"]["useafter"] == "<SPECIFY>"
    # No validation errors against the simple schema:
    with asdf.AsdfFile(asdf_file.tree) as af:
        _validate_asdf(af, "http://stsci.edu/schemas/stpipe/step_config-1.0.0")

    # Expect failures from the metadata schema while <SPECIFY>
    # is still present:
    with pytest.raises(asdf.ValidationError):
        with asdf.AsdfFile(asdf_file.tree) as af:
            _validate_asdf(af, "http://stsci.edu/schemas/stpipe/step_config_with_metadata-1.0.0")

    asdf_file["meta"]["author"] = "Yours Truly"
    asdf_file["meta"]["instrument"]["name"] = "EYEBALLS"
    asdf_file["meta"]["origin"] = "The depths"
    asdf_file["meta"]["pedigree"] = "GROUND"
    asdf_file["meta"]["reftype"] = "pars-pipelineclass"
    asdf_file["meta"]["telescope"] = "BINOCULARS"
    asdf_file["meta"]["useafter"] = "2020-11-20T00:00:00"
    # No validation errors now:
    with asdf.AsdfFile(asdf_file.tree) as af:
        _validate_asdf(af, "http://stsci.edu/schemas/stpipe/step_config_with_metadata-1.0.0")


def test_step_config_from_asdf(config):
    asdf_file = config.to_asdf()
    assert StepConfig.from_asdf(asdf_file) == config

    with pytest.raises(asdf.ValidationError):
        StepConfig.from_asdf(asdf.AsdfFile())


def test_step_config_from_legacy_asdf(config):
    parameters = {
        "class": config.class_name,
        "name": config.name,
    }
    parameters.update(config.parameters)

    steps = {}
    for step in config.steps:
        # The older config files don't include "name" or "class"
        # for the child steps.
        steps[step.name] = step.parameters

    parameters["steps"] = steps

    asdf_file = asdf.AsdfFile({"parameters": parameters}, custom_schema="http://stsci.edu/schemas/stpipe/step_config-0.1.0")

    # We lose the step class names so the configs won't be equal:
    legacy_config = StepConfig.from_asdf(asdf_file)
    assert legacy_config.name == config.name
    assert legacy_config.class_name == config.class_name
    assert legacy_config.parameters == config.parameters
    assert len(legacy_config.steps) == len(config.steps)
    for legacy_step, step in zip(legacy_config.steps, config.steps):
        assert legacy_step.name == step.name
        assert legacy_step.class_name is None
        assert legacy_step.parameters == step.parameters
        assert legacy_step.steps == step.steps


def test_export_config_step():
    step = WithDefaultsStep(name="Ronald", par1="override par1 value")
    config = export_config(step)

    assert config.class_name == "jwst.stpipe.tests.steps.WithDefaultsStep"
    assert config.name == "Ronald"

    assert config.parameters["par1"] == "override par1 value"
    assert config.parameters["par2"] == "default par2 value"
    assert config.parameters["par3"] == "default par3 value"
    assert config.parameters["par4"] == "default par4 value"

    # Not going to check all superclass parameters, just one to make
    # sure that Step was consulted:
    assert config.parameters["input_dir"] == ""

    assert len(config.steps) == 0


def test_export_config_pipeline():
    pipeline = MakeListPipeline(
        name="Dorothy",
        par1="override the atomizer",
        steps={"make_list": {"par1": 3.14159, "par2": "frabjous"}},
    )
    config = export_config(pipeline)

    assert config.class_name == "jwst.stpipe.tests.steps.MakeListPipeline"
    assert config.name == "Dorothy"

    assert config.parameters["par1"] == "override the atomizer"

    # Not going to check all superclass parameters, just one to make
    # sure that Step was consulted:
    assert config.parameters["input_dir"] == ""

    assert len(config.steps) == 1
    step = config.steps[0]

    assert step.class_name == "jwst.stpipe.tests.steps.MakeListStep"
    assert step.name == "make_list"

    assert step.parameters["par1"] == 3.14159
    assert step.parameters["par2"] == "frabjous"
    assert step.parameters["par3"] is False

    # Check a base class parameter on the step too:
    assert step.parameters["input_dir"] == ""
