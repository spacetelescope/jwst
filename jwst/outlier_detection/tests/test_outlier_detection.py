import numpy as np
import pytest
from gwcs.wcs import WCS
from numpy.testing import assert_allclose, assert_array_equal
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.tests import helpers
from jwst.outlier_detection.utils import (
    _flag_resampled_model_crs,
    median_with_resampling,
    median_without_resampling,
)


def test_flag_cr(sci_blot_image_pair):
    """Test the flag_cr function.  Test logic, not the actual noise model."""
    sci, blot = sci_blot_image_pair
    assert_array_equal(sci.dq, 0)

    # Drop some CRs on the science array
    sci.data[3, 3] += 100
    sci.data[3, 7] += 1e3
    sci.data[7, 3] += 1e4
    sci.data[7, 7] += 1e5

    # run flag_cr() which updates in-place.  Copy sci first.
    data_copy = sci.data.copy()
    _flag_resampled_model_crs(
        sci,
        blot.data,
        None,
        5.0,
        4.0,
        1.2,
        0.7,
        0,
    )

    # Make sure science data array is unchanged after flag_cr()
    # except outliers are NaN
    dnu = (sci.dq & helpers.OUTLIER_DO_NOT_USE).astype(bool)
    assert np.all(np.isnan(sci.data[dnu]))
    assert_allclose(sci.data[~dnu], data_copy[~dnu])

    # Verify that both DQ flags are set in the DQ array for all outliers
    assert sci.dq[3, 3] == helpers.OUTLIER_DO_NOT_USE
    assert sci.dq[3, 7] == helpers.OUTLIER_DO_NOT_USE
    assert sci.dq[7, 3] == helpers.OUTLIER_DO_NOT_USE
    assert sci.dq[7, 7] == helpers.OUTLIER_DO_NOT_USE

    # Verify the source wasn't flagged
    assert sci.dq[10, 10] == datamodels.dqflags.pixel["GOOD"]


def test_same_median_on_disk(three_sci_as_asn, tmp_cwd):
    """Test creation of median on disk vs in memory"""
    lib_on_disk = ModelLibrary(three_sci_as_asn, on_disk=True)
    lib_in_memory = ModelLibrary(three_sci_as_asn, on_disk=False)

    # make this test meaningful w.r.t. handling of weights
    with lib_on_disk, lib_in_memory:
        for lib in [lib_on_disk, lib_in_memory]:
            for model in lib:
                model.var_rnoise = np.ones_like(model.data)
                model.var_rnoise[4, 9] = 2.0
                lib.shelve(model, modify=True)

    # 32-bit floats are 4 bytes each, min buffer size is one row of 20 pixels
    # arbitrarily use 5 times that
    buffer_size = 4 * 20 * 5
    median_on_disk, _ = median_without_resampling(
        lib_on_disk,
        0.7,
        "ivm",
        "~DO_NOT_USE",
        buffer_size=buffer_size,
    )
    median_in_memory, _ = median_without_resampling(
        lib_in_memory,
        0.7,
        "ivm",
        "~DO_NOT_USE",
        buffer_size=buffer_size,
    )

    # Make sure the high-variance (low-weight) pixel is set to NaN
    assert np.isnan(median_in_memory[4, 9])

    # Make sure the median library is the same for on-disk and in-memory
    assert np.allclose(median_on_disk, median_in_memory, equal_nan=True)


def test_drizzle_and_median_with_resample(three_sci_as_asn, tmp_cwd):
    lib = ModelLibrary(three_sci_as_asn, on_disk=False)

    resamp = helpers.make_resamp(lib)
    median, wcs = median_with_resampling(lib, resamp, 0.7)

    assert isinstance(wcs, WCS)
    assert median.shape == (34, 34)

    resamp.single = False
    with pytest.raises(ValueError):
        # ensure failure if try to call when resamp.single is False
        median_with_resampling(lib, resamp, 0.7, save_intermediate_results=True)


@pytest.mark.parametrize("mode", [None, "unknown"])
def test_guess_mode_assigned(caplog, mode):
    input_model = datamodels.ImageModel()

    # Assign a "None" mode -
    # processing should gracefully skip
    step = OutlierDetectionStep()
    step.mode = mode
    result = step.run(input_model)

    assert result.meta.cal_step.outlier_detection == "SKIPPED"
    assert isinstance(result, datamodels.ImageModel)

    # If mode was unrecognized, error message is issued
    if mode is None:
        assert "ERROR" not in caplog.text
    else:
        assert "ERROR" in caplog.text

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.outlier_detection is None


def test_skip_unknown_mode_file(tmp_path, caplog):
    input_model = datamodels.ImageModel()
    input_file = str(tmp_path / "test.fits")
    input_model.save(input_file)
    result = OutlierDetectionStep.call(input_file)

    # Step is skipped with an error message
    assert result.meta.cal_step.outlier_detection == "SKIPPED"
    assert isinstance(result, datamodels.ImageModel)
    assert "ERROR" in caplog.text

    # Input is not modified
    with datamodels.open(input_file) as input_model:
        assert result is not input_model
        assert input_model.meta.cal_step.outlier_detection is None


def test_skip_unknown_mode_image_model():
    input_model = datamodels.ImageModel()
    result = OutlierDetectionStep.call(input_model)

    # Step is skipped
    assert result.meta.cal_step.outlier_detection == "SKIPPED"

    # Input is not modified
    assert result is not input_model
    assert input_model.meta.cal_step.outlier_detection is None


def test_skip_unknown_mode_container():
    input_model = datamodels.ImageModel()
    container = ModelContainer([input_model])
    result = OutlierDetectionStep.call(container)

    # Step is skipped
    assert result[0].meta.cal_step.outlier_detection == "SKIPPED"

    # Input is not modified
    assert result[0] is not input_model
    assert input_model.meta.cal_step.outlier_detection is None
