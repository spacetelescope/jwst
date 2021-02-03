import glob
import os

import asdf
import yaml

from jwst.stpipe.integration import get_steps, SCHEMAS_PATH
from jwst.stpipe import Step, Pipeline
from jwst.stpipe.utilities import import_class

import jwst.pipeline
import jwst.step


def test_get_steps():
    tuples = get_steps()

    assert {t[0].split(".")[-1] for t in tuples} == set(jwst.step.__all__ + jwst.pipeline.__all__)

    for class_name, class_alias, is_pipeline in tuples:
        step_class = import_class(class_name)
        assert issubclass(step_class, Step)
        assert step_class.class_alias == class_alias
        if is_pipeline:
            assert issubclass(step_class, Pipeline)


def test_asdf_extension():
    for schema_path in glob.glob(os.path.join(SCHEMAS_PATH, "**/*.yaml"), recursive=True):
        with open(schema_path) as f:
            yaml_schema = yaml.safe_load(f.read())
            asdf_schema = asdf.schema.load_schema(yaml_schema["id"])
            assert asdf_schema["id"] == yaml_schema["id"]
