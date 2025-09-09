from stdatamodels.jwst import datamodels

from jwst.srctype.srctype_step import SourceTypeStep


def test_step_complete():
    model = datamodels.ImageModel()
    model.meta.exposure.type = "NRS_IFU"
    model.meta.target.source_type = "UNKNOWN"

    result = SourceTypeStep.call(model)

    # Step is complete
    assert result.meta.cal_step.srctype == "COMPLETE"
    assert result.meta.target.source_type == "EXTENDED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.srctype is None
    assert model.meta.target.source_type == "UNKNOWN"


def test_srctype_provided():
    model = datamodels.ImageModel()
    model.meta.exposure.type = "NRS_IFU"
    model.meta.target.source_type = "UNKNOWN"

    result = SourceTypeStep.call(model, source_type="POINT")

    # Step is complete
    assert result.meta.cal_step.srctype == "COMPLETE"
    assert result.meta.target.source_type == "POINT"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.srctype is None
    assert model.meta.target.source_type == "UNKNOWN"
