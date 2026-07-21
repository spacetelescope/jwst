from pathlib import Path

import numpy as np

from jwst.pipeline import DarkPipeline
from jwst.pipeline.tests.helpers import make_miri_ramp_model


def test_dark_pipeline_model(tmp_path):
    input_model = make_miri_ramp_model()
    input_model_copy = input_model.copy()
    DarkPipeline.call(input_model, save_results=True, output_dir=str(tmp_path))

    # Check for expected output
    assert (tmp_path / "test_miri_dark.fits").exists()

    # Input is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step.instance == input_model_copy.meta.cal_step.instance

    input_model.close()
    input_model_copy.close()


def test_dark_pipeline_asn(tmp_cwd):
    prodname = "test_dark"
    modname = f"{prodname}_target_1.fits"
    input_model = make_miri_ramp_model()
    input_model.save(modname)
    asn = {
        "asn_id": "c1000",
        "target": "t001",
        "asn_pool": f"{prodname}_pool.csv",
        "products": [{"name": prodname, "members": [{"expname": modname, "exptype": "science"}]}],
    }
    DarkPipeline.call(asn, save_results=True, output_file="test_miri")

    # Check for expected output
    assert Path("test_miri_dark.fits").exists()

    input_model.close()
