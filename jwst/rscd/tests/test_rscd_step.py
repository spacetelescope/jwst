from jwst.rscd.rscd_step import RscdStep


def test_step_complete(create_miri_model):
    model = create_miri_model()
    result = RscdStep.call(model)

    # Step is marked complete
    assert result.meta.cal_step.rscd == "COMPLETE"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.rscd is None


def test_skip_nonmiri(caplog, create_miri_model):
    model = create_miri_model()
    model.meta.instrument.detector = "NRS1"
    result = RscdStep.call(model)

    # Step is skipped
    assert "only for MIRI data" in caplog.text
    assert result.meta.cal_step.rscd == "SKIPPED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.rscd is None


def test_skip_missing_reffile(caplog, create_miri_model):
    model = create_miri_model()
    result = RscdStep.call(model, override_rscd="N/A")

    # Step is skipped
    assert "No RSCD reference file found" in caplog.text
    assert result.meta.cal_step.rscd == "SKIPPED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.rscd is None


def test_skip_empty_param(caplog, create_miri_model):
    # make a model with an unrecognized subarray
    model = create_miri_model()
    model.meta.subarray.name = "FULLP"
    result = RscdStep.call(model)

    # Step is skipped
    assert "READPATT, SUBARRAY combination not found" in caplog.text
    assert result.meta.cal_step.rscd == "SKIPPED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.rscd is None
