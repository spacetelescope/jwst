"""Test utility funcs"""
from stpipe.utilities import resolve_step_class_alias

from jwst.stpipe.utilities import all_steps
import jwst.pipeline
import jwst.step


# All steps available to users should be represented in one
# of these two modules:
KNOWN_STEPS = set(jwst.pipeline.__all__ + jwst.step.__all__)


def test_all_steps():
    """Test finding all defined steps and pipelines"""
    found_steps = all_steps()
    diff = KNOWN_STEPS.symmetric_difference(found_steps)
    assert not diff, f'Steps not accounted for. Confirm and check suffix and CRDS calpars.\n{diff}'


def test_resolve_step_class_alias():
    for class_name in jwst.pipeline.__all__:
        pipeline_class = getattr(jwst.pipeline, class_name)
        full_class_name = f"jwst.pipeline.{class_name}"
        if pipeline_class.class_alias is not None:
            assert resolve_step_class_alias(pipeline_class.class_alias) == full_class_name
        assert resolve_step_class_alias(full_class_name) == full_class_name
