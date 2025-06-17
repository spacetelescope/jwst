import pytest
from jwst.outlier_detection._fileio import _save_intermediate_output
from jwst.datamodels import ImageModel  # type: ignore[attr-defined]
from jwst.step import OutlierDetectionStep
import os
import numpy as np
from functools import partial


SLIT_ID = "slit007"
ASN_ID = "a008"
INPUT_FILENAME = "foo_outlier_s2d.fits"
SUFFIX = "test"


@pytest.fixture(scope="module")
def model():
    model = ImageModel(data=np.zeros((10, 10)))
    model.meta.filename = INPUT_FILENAME
    return model


@pytest.fixture(scope="module")
def make_output_path():
    return OutlierDetectionStep().make_output_path


@pytest.mark.parametrize("asn_id", [None, ASN_ID])
@pytest.mark.parametrize("slit_id", [None, SLIT_ID])
def test_save(tmp_cwd, model, make_output_path, asn_id, slit_id):
    this_model = model.copy()
    if slit_id is not None:
        this_model.name = slit_id
    make_output_path = partial(make_output_path, asn_id=asn_id)
    _save_intermediate_output(this_model, SUFFIX, make_output_path)

    stem = model.meta.filename.split("_")[0]
    inputs = [val for val in [stem, asn_id, slit_id, SUFFIX] if val is not None]
    expected_filename = "_".join(inputs) + ".fits"
    assert os.path.isfile(expected_filename)
