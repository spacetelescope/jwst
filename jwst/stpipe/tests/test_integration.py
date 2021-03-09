from stpipe.utilities import import_class

from jwst.stpipe.integration import get_steps
from jwst.stpipe import Step, Pipeline

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
