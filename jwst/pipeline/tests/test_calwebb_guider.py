import numpy as np

from jwst.guider_cds.tests.helpers import make_guider_image
from jwst.pipeline import GuiderPipeline


def test_guider_pipeline(tmp_path):
    input_model = make_guider_image()
    input_model_copy = input_model.copy()
    GuiderPipeline.call(input_model, save_results=True, output_dir=str(tmp_path))

    # Check for expected output
    assert (tmp_path / "test_guider_cal.fits").exists()

    # Input model is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step._instance == input_model_copy.meta.cal_step._instance

    input_model.close()
    input_model_copy.close()
