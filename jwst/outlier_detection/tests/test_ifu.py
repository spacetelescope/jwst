import numpy as np
import pytest

from jwst.datamodels import ModelContainer
from jwst.outlier_detection.outlier_detection_step import OutlierDetectionStep


def test_ifu_one_exposure(miri_ifu_rate):
    input_model = miri_ifu_rate
    result = OutlierDetectionStep.call([input_model])

    # Step is skipped
    assert result[0].meta.cal_step.outlier_detection == "SKIPPED"

    # Input is not modified
    assert result[0] is not input_model
    assert input_model.meta.cal_step.outlier_detection is None


@pytest.mark.parametrize("dataset,expected", [("miri_ifu_rate", 20), ("nirspec_ifu_rate", 23)])
def test_ifu_with_outliers(caplog, request, dataset, expected):
    input_model_1 = request.getfixturevalue(dataset)
    input_model_2 = input_model_1.copy()
    input_models = ModelContainer([input_model_1, input_model_2])
    result = OutlierDetectionStep.call(input_models, kernel_size="7 7")

    # Step is complete
    assert result[0].meta.cal_step.outlier_detection == "COMPLETE"

    # Both images flag the same number of outliers
    assert caplog.text.count(f"Total # pixels flagged as outliers: {expected}") == 2
    for model in result:
        assert np.sum(np.isnan(model.data)) == expected
        assert np.sum(np.isnan(model.err)) == expected
        assert np.sum(model.dq & 1 > 0) == expected

    # Input is not modified
    for model in input_models:
        assert result[0] is not model
        assert result[1] is not model
        assert model.meta.cal_step.outlier_detection is None
        assert np.sum(np.isnan(model.data)) == 0
        assert np.sum(np.isnan(model.err)) == 0
        assert np.sum(model.dq & 1 > 0) == 0


def test_ifu_second_check(caplog, miri_ifu_rate):
    input_model_1 = miri_ifu_rate
    input_model_2 = miri_ifu_rate.copy()
    result = OutlierDetectionStep.call(
        [input_model_1, input_model_2], kernel_size="7 7", ifu_second_check=True
    )

    # Step is complete
    assert result[0].meta.cal_step.outlier_detection == "COMPLETE"
    assert caplog.text.count("Total # pixels flagged as outliers: 20") == 2

    # No additional outliers in second check
    assert "second check: 0" in caplog.text

    # Input is not modified
    assert result[0] is not input_model_1
    assert result[1] is not input_model_2
    assert input_model_1.meta.cal_step.outlier_detection is None
    assert input_model_2.meta.cal_step.outlier_detection is None


def test_ifu_even_kernel(caplog, miri_ifu_rate):
    input_model_1 = miri_ifu_rate
    input_model_2 = miri_ifu_rate.copy()
    OutlierDetectionStep.call([input_model_1, input_model_2], kernel_size="6 6")

    # kernel is sized up to "7 7"
    assert caplog.text.count("Increasing number by 1") == 2

    # Both images still flag 20 outliers
    assert caplog.text.count("Total # pixels flagged as outliers: 20") == 2


def test_ifu_save_intermediate(tmp_path, miri_ifu_rate):
    input_model_1 = miri_ifu_rate
    input_model_1.meta.filename = "test1.fits"
    input_model_2 = miri_ifu_rate.copy()
    input_model_2.meta.filename = "test2.fits"
    OutlierDetectionStep.call(
        [input_model_1, input_model_2],
        suffix="outlier_detection",
        output_dir=str(tmp_path),
        save_results=True,
        save_intermediate_results=True,
    )
    assert (tmp_path / "test1_outlier_detection.fits").exists()
    assert (tmp_path / "test2_outlier_detection.fits").exists()
    assert (tmp_path / "test1_mirifushort_outlier_output.fits").exists()
