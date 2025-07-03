"""Test utility funcs"""

from stpipe.utilities import resolve_step_class_alias
from jwst.stpipe import record_step_status, query_step_status

from jwst.stpipe.utilities import all_steps, NOT_SET
import jwst.pipeline
import jwst.step
from stdatamodels.jwst import datamodels
from jwst.stpipe.utilities import invariant_filename

from jwst import datamodels as dm


# All steps available to users should be represented in one
# of these two modules:
KNOWN_STEPS = set(jwst.pipeline.__all__ + jwst.step.__all__)


def test_all_steps():
    """Test finding all defined steps and pipelines"""
    found_steps = all_steps()
    diff = KNOWN_STEPS.symmetric_difference(found_steps)
    assert not diff, f"Steps not accounted for. Confirm and check suffix and CRDS calpars.\n{diff}"


def test_resolve_step_class_alias():
    for class_name in jwst.pipeline.__all__:
        pipeline_class = getattr(jwst.pipeline, class_name)
        full_class_name = f"jwst.pipeline.{class_name}"
        if pipeline_class.class_alias is not None:
            assert resolve_step_class_alias(pipeline_class.class_alias) == full_class_name
        assert resolve_step_class_alias(full_class_name) == full_class_name


def test_record_query_step_status():
    # single model case
    model = dm.MultiSpecModel()
    record_step_status(model, "test_step", success=True)
    assert query_step_status(model, "test_step") == "COMPLETE"

    # modelcontainer case where all skipped
    model2 = dm.ModelContainer()
    model2.append(dm.MultiSpecModel())
    model2.append(dm.MultiSpecModel())
    record_step_status(model2, "test_step", success=False)
    assert model2[0].meta.cal_step._instance["test_step"] == "SKIPPED"
    assert model2[1].meta.cal_step._instance["test_step"] == "SKIPPED"
    assert query_step_status(model2, "test_step") == "SKIPPED"

    # modelcontainer case with at least one complete
    # currently the zeroth is checked, so this should return SKIPPED
    model2.append(dm.MultiSpecModel())
    model2[2].meta.cal_step._instance["test_step"] = "COMPLETE"
    assert query_step_status(model2, "test_step") == "SKIPPED"

    # test query not set
    model3 = dm.MultiSpecModel()
    assert query_step_status(model3, "test_step") == NOT_SET


def change_name_func(model):
    model.meta.filename = "changed"
    model.meta.cal_step.pixel_replace = "COMPLETE"
    return model


def test_invariant_filename():
    # Make sure the change_name_func changes the name and has side effects
    # (here, setting a status variable, but normally, actually saving the file)
    input_model = datamodels.IFUImageModel()
    input_model.meta.filename = "test1.fits"
    change_name_func(input_model)
    assert input_model.meta.filename == "changed"
    assert input_model.meta.cal_step.pixel_replace == "COMPLETE"

    # When the function is wrapped with invariant_filename,
    # the filename is not changed, but the side effect still happens
    input_model = datamodels.IFUImageModel()
    input_model.meta.filename = "test2.fits"
    invariant_save_func = invariant_filename(change_name_func)
    output_model = invariant_save_func(input_model)
    assert output_model.meta.filename == "test2.fits"
    assert output_model.meta.cal_step.pixel_replace == "COMPLETE"

    # The output model is not a copy - the name is reset in place
    assert output_model is input_model
