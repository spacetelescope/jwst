import numpy as np

from jwst.pipeline import DarkPipeline
from jwst.pipeline.tests.helpers import make_miri_ramp_model


def test_dark_pipeline(tmp_path):
    input_model = make_miri_ramp_model()
    input_model_copy = input_model.copy()
    DarkPipeline.call(input_model, save_results=True, output_dir=str(tmp_path))

    # Check for expected output
    assert (tmp_path / "test_miri_dark.fits").exists()

    # Input is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step._instance == input_model_copy.meta.cal_step._instance

    input_model.close()
    input_model_copy.close()
