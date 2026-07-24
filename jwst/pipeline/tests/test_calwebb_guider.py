from pathlib import Path

import numpy as np

from jwst.guider_cds.tests.helpers import make_guider_image
from jwst.pipeline import GuiderPipeline


def test_guider_pipeline_model(tmp_path):
    input_model = make_guider_image()
    input_model_copy = input_model.copy()
    GuiderPipeline.call(input_model, save_results=True, output_dir=str(tmp_path))

    # Check for expected output
    assert (tmp_path / "test_guider_cal.fits").exists()

    # Input model is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step.instance == input_model_copy.meta.cal_step.instance

    input_model.close()
    input_model_copy.close()


def test_guider_pipeline_asn(tmp_cwd):
    prodname = "test_guider"
    modname = f"{prodname}_target_1.fits"
    input_model = make_guider_image()
    input_model.save(modname)
    asn = {
        "asn_id": "c1000",
        "target": "t001",
        "asn_pool": f"{prodname}_pool.csv",
        "products": [{"name": prodname, "members": [{"expname": modname, "exptype": "science"}]}],
    }
    GuiderPipeline.call(asn, save_results=True, output_file="test_guider")

    # Check for expected output
    assert Path("test_guider_cal.fits").exists()

    input_model.close()
