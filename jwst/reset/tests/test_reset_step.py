from jwst.reset.reset_step import ResetStep


def test_step_complete(make_rampmodel):
    model = make_rampmodel()
    result = ResetStep.call(model)

    # Step is complete
    assert result.meta.cal_step.reset == "COMPLETE"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.reset is None


def test_skip_nonmiri(caplog, make_rampmodel):
    model = make_rampmodel()
    model.meta.instrument.detector = "NRS1"
    result = ResetStep.call(model)

    # Step is skipped
    assert "only for MIRI data" in caplog.text
    assert result.meta.cal_step.reset == "SKIPPED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.reset is None


def test_skip_missing_reffile(caplog, make_rampmodel):
    model = make_rampmodel()
    result = ResetStep.call(model, override_reset="N/A")

    # Step is skipped
    assert "No RESET reference file found" in caplog.text
    assert result.meta.cal_step.reset == "SKIPPED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.reset is None
