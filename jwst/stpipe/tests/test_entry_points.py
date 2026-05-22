import jwst
import jwst.pipeline
import jwst.step
from jwst.stpipe import Pipeline, Step
from jwst.stpipe.integration import get_steps


def test_get_steps():
    step_info_by_class = {s[0]: s for s in get_steps()}
    class_paths = set()

    for class_name in jwst.pipeline.__all__:
        pipeline_class = getattr(jwst.pipeline, class_name)
        class_path = f"jwst.pipeline.{class_name}"
        _, class_alias, is_pipeline = step_info_by_class[class_path]

        assert class_alias == pipeline_class.class_alias
        assert is_pipeline is True
        assert issubclass(pipeline_class, Step)
        assert issubclass(pipeline_class, Pipeline)
        class_paths.add(class_path)

    for class_name in jwst.step.__all__:
        step_class = getattr(jwst.step, class_name)
        class_path = f"jwst.step.{class_name}"
        _, class_alias, is_pipeline = step_info_by_class[class_path]
        assert class_alias == step_class.class_alias
        assert is_pipeline is False
        assert issubclass(pipeline_class, Step)
        class_paths.add(class_path)

    assert class_paths == set(step_info_by_class.keys())
