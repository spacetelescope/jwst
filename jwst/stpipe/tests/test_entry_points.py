import jwst
import jwst.pipeline
import jwst.step
from jwst.stpipe.integration import get_steps


def test_get_steps():
    step_info_by_class = {s[0]: s for s in get_steps()}

    for class_name in jwst.pipeline.__all__:
        pipeline_class = getattr(jwst.pipeline, class_name)
        _, class_alias, is_pipeline = step_info_by_class[f"jwst.pipeline.{class_name}"]

        assert class_alias == pipeline_class.class_alias
        assert is_pipeline is True

    for class_name in jwst.step.__all__:
        step_class = getattr(jwst.step, class_name)
        _, class_alias, is_pipeline = step_info_by_class[f"jwst.step.{class_name}"]
        assert class_alias == step_class.class_alias
        assert is_pipeline is False
