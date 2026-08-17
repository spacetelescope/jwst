import numpy as np
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline


def test_load_single_datamodel_as_container():
    model = datamodels.ImageModel(data=np.ones((2, 2)))
    model.meta.filename = "jw00001_cal.fits"
    pipeline = Spec3Pipeline()
    result = pipeline._load_input_as_container(model, ["science", "background"])
    try:
        assert isinstance(result, ModelContainer)
        assert len(result) == 1
        assert result.asn_table["products"][0]["name"] == "jw00001"
    finally:
        result.close()
        model.close()


def test_load_single_fits_as_container(tmp_path):
    path = tmp_path / "jw00002_cal.fits"
    with datamodels.ImageModel(data=np.ones((2, 2))) as model:
        model.save(path)
    pipeline = Spec3Pipeline()
    result = pipeline._load_input_as_container(str(path), ["science", "background"])
    try:
        assert isinstance(result, ModelContainer)
        assert len(result) == 1
        assert result.asn_table["products"][0]["name"] == "jw00002"
    finally:
        result.close()
