from stdatamodels.jwst import datamodels

from jwst.extract_2d.extract_2d_step import Extract2dStep


def test_skip_nrs_mirror(caplog):
    model = datamodels.ImageModel()
    model.meta.exposure.type = "NRS_LAMP"
    model.meta.instrument.grating = "MIRROR"

    # Step is skipped
    result = Extract2dStep.call(model)
    assert result.meta.cal_step.extract_2d == "SKIPPED"
    assert "grating=MIRROR not supported" in caplog.text

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.extract_2d is None


def test_skip_unsupported(caplog):
    model = datamodels.ImageModel()
    model.meta.exposure.type = "NRC_IMAGE"

    # Step is skipped
    result = Extract2dStep.call(model)
    assert result.meta.cal_step.extract_2d == "SKIPPED"
    assert "EXP_TYPE NRC_IMAGE not supported" in caplog.text

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.extract_2d is None
