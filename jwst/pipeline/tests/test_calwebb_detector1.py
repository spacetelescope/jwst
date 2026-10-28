import numpy as np
from stdatamodels.jwst import datamodels

from jwst.pipeline import Detector1Pipeline
from jwst.pipeline.tests.helpers import make_miri_ramp_model


def make_miri_dark(tmp_path, input_model):
    dark_model = input_model.copy()
    dark_model.meta.description = "lorem ipsum"
    dark_model.meta.reftype = "DARK"
    dark_model.meta.author = "Dolor"
    dark_model.meta.pedigree = "unit test"
    dark_model.meta.useafter = "2022-04-01T00:00:00"
    dark_model.meta.telescope = "JWST"
    dark_model = datamodels.DarkMIRIModel(dark_model)
    dark_model.save(tmp_path / "test_miri_dark.fits")


def test_detector1_miri(tmp_path):
    input_model = make_miri_ramp_model()
    input_model_copy = input_model.copy()
    make_miri_dark(tmp_path, input_model)

    Detector1Pipeline.call(
        input_model,
        save_results=True,
        output_dir=str(tmp_path),
        steps={"dark_current": {"override_dark": "test_miri_dark.fits"}},
    )

    # Check for expected output
    assert (tmp_path / "test_miri_rate.fits").exists()
    assert (tmp_path / "test_miri_rateints.fits").exists()

    # Input is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step.instance == input_model_copy.meta.cal_step.instance

    input_model.close()
    input_model_copy.close()
