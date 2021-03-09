from stpipe.entry_points import get_steps

import jwst.pipeline
import jwst.step
import jwst


def test_get_steps():
    step_info_by_class = {s.class_name: s for s in get_steps()}

    for class_name in jwst.pipeline.__all__:
        pipeline_class = getattr(jwst.pipeline, class_name)
        step_info = step_info_by_class[f"jwst.pipeline.{class_name}"]
        assert step_info.class_alias == pipeline_class.class_alias
        assert step_info.is_pipeline is True
        assert step_info.package_name == "jwst"
        assert step_info.package_version == jwst.__version__

    for class_name in jwst.step.__all__:
        step_class = getattr(jwst.step, class_name)
        step_info = step_info_by_class[f"jwst.step.{class_name}"]
        assert step_info.class_alias == step_class.class_alias
        assert step_info.is_pipeline is False
        assert step_info.package_name == "jwst"
        assert step_info.package_version == jwst.__version__
