"""Test utility funcs"""

import pytest
from stdatamodels.jwst import datamodels

import jwst.pipeline
import jwst.step
from jwst import datamodels as dm
from jwst.stpipe import query_step_status, record_step_status
from jwst.stpipe.utilities import NOT_SET, all_steps, invariant_filename, summary_step_status

# All steps available to users should be represented in one
# of these two modules:
KNOWN_STEPS = set(jwst.pipeline.__all__ + jwst.step.__all__)


def test_all_steps():
    """Test finding all defined steps and pipelines"""
    found_steps = all_steps()
    diff = KNOWN_STEPS.symmetric_difference(found_steps)
    assert not diff, f"Steps not accounted for. Confirm and check suffix and CRDS calpars.\n{diff}"


@pytest.mark.parametrize("fail_status", [None, "SKIPPED", "FAILED"])
def test_record_query_step_status(fail_status):
    # single model case
    model = dm.MultiSpecModel()
    record_step_status(model, "test_step", success=True)
    assert query_step_status(model, "test_step") == "COMPLETE"

    # modelcontainer case where all failed
    model2 = dm.ModelContainer()
    model2.append(dm.MultiSpecModel())
    model2.append(dm.MultiSpecModel())
    record_step_status(model2, "test_step", success=False, status=fail_status)
    expected = fail_status if fail_status is not None else "FAILED"
    assert model2[0].meta.cal_step.instance["test_step"] == expected
    assert model2[1].meta.cal_step.instance["test_step"] == expected
    assert query_step_status(model2, "test_step") == expected

    # modelcontainer case with at least one complete
    # currently the zeroth is checked, so this should return FAILED
    model2.append(dm.MultiSpecModel())
    model2[2].meta.cal_step.instance["test_step"] = "COMPLETE"
    assert query_step_status(model2, "test_step") == expected

    # test query not set
    model3 = dm.MultiSpecModel()
    assert query_step_status(model3, "test_step") == NOT_SET

    # test library
    model4 = dm.ModelLibrary([model])
    assert query_step_status(model4, "test_step") == "COMPLETE"
    model5 = dm.ModelLibrary(model2)
    assert query_step_status(model5, "test_step") == expected
    model6 = dm.ModelLibrary([model3])
    assert query_step_status(model6, "test_step") == NOT_SET


def test_record_status_invalid_value():
    model = dm.MultiSpecModel()
    with pytest.raises(ValueError, match="BAD_VALUE not in allowed values"):
        record_step_status(model, "test_step", status="BAD_VALUE")


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


@pytest.mark.parametrize(
    "status_values,expected",
    [
        (("COMPLETE", "COMPLETE"), "COMPLETE"),
        (("COMPLETE", "SKIPPED"), "COMPLETE"),
        (("COMPLETE", "FAILED"), "COMPLETE"),
        (("COMPLETE", "BAD"), "COMPLETE"),
        (("FAILED", "FAILED"), "FAILED"),
        (("FAILED", "SKIPPED"), "FAILED"),
        (("FAILED", "BAD"), "FAILED"),
        (("SKIPPED", "SKIPPED"), "SKIPPED"),
        (("SKIPPED", "BAD"), "SKIPPED"),
    ],
)
def test_summary_step_status(status_values, expected):
    assert summary_step_status(status_values) == expected
