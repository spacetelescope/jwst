from stpipe.utilities import import_class

from jwst.stpipe.integration import get_steps
from jwst.stpipe import Step, Pipeline
from jwst.stpipe.utilities import NON_STEPS

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


def _get_subclass_names(class_, ignore=None, ignore_tests=True):
    """
    Iterate through names of all subclasses (recursively) of a class ignoring
    all subclasses with names in ignore (and ignoring any subclasses of
    ignored classes). Also optionally ignore all classes defined in tests modules.
    """
    if ignore is None:
        ignore = set()
    for subclass in class_.__subclasses__():
        if subclass.__name__ in ignore:
            continue
        if ignore_tests and (
            subclass.__module__ == "local" or "tests" in subclass.__module__.split(".")
        ):
            continue
        yield subclass.__name__
        yield from _get_subclass_names(subclass, ignore, ignore_tests)


def test_all_steps_in_step_all():
    step_names = set(_get_subclass_names(Step, set(["JwstPipeline"])))
    # MasterBackgroundMosStep is a pipeline that is treated like a step
    step_names.add("MasterBackgroundMosStep")
    # also ignore any 'non-steps'
    step_names -= set(NON_STEPS)
    assert set(jwst.step.__all__) == step_names


def test_all_pipelines_in_pipeline_all():
    pipeline_names = set(_get_subclass_names(Pipeline))
    # MasterBackgroundMosStep is a pipeline that is treated like a step
    pipeline_names.remove("MasterBackgroundMosStep")
    assert set(jwst.pipeline.__all__) == pipeline_names
