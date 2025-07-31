import os
from glob import glob

import numpy as np
import pytest
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.tests import helpers
from jwst.stpipe.utilities import query_step_status


@pytest.mark.parametrize("do_resample", [True, False])
def test_outlier_step_no_outliers(mirimage_three_sci, do_resample, tmp_cwd):
    """Test whole step, no outliers"""
    container = ModelContainer(list(mirimage_three_sci))
    container[0].var_rnoise[10, 10] = 1e9
    pristine = ModelContainer([m.copy() for m in container])
    OutlierDetectionStep.call(container, in_memory=True, resample_data=do_resample)

    # Make sure nothing changed in SCI and DQ arrays
    for image, uncorrected in zip(pristine, container):
        np.testing.assert_allclose(image.data, uncorrected.data)
        np.testing.assert_allclose(image.dq, uncorrected.dq)


def test_outlier_step_weak_cr_imaging(mirimage_three_sci, tmp_cwd):
    """Test whole step with an outlier including saving intermediate and results files"""
    container = ModelLibrary(mirimage_three_sci)

    # Drop a CR on the science array
    with container:
        zeroth = container.borrow(0)
        zeroth.data[12, 12] += 1
        container.shelve(zeroth)

    # Verify that intermediate files are removed
    OutlierDetectionStep.call(container)
    i2d_files = glob(os.path.join(tmp_cwd, "*i2d.fits"))
    median_files = glob(os.path.join(tmp_cwd, "*median.fits"))
    assert len(i2d_files) == 0
    assert len(median_files) == 0

    # Save all the data into a separate array before passing into step
    data_as_cube = list(
        container.map_function(lambda model, index: model.data.copy(), modify=False)
    )

    result = OutlierDetectionStep.call(container, save_results=True, save_intermediate_results=True)

    with result:
        for i, r in enumerate(result):
            # Make sure nothing changed in SCI array except outliers are NaN
            dnu = (r.dq & helpers.OUTLIER_DO_NOT_USE).astype(bool)
            assert np.all(np.isnan(r.data[dnu]))
            assert_allclose(data_as_cube[i][~dnu], r.data[~dnu])

            # Verify source is not flagged
            assert r.dq[7, 7] == datamodels.dqflags.pixel["GOOD"]
            result.shelve(r, modify=False)

    # Verify CR is flagged
    with result:
        zeroth = result.borrow(0)
        assert zeroth.dq[12, 12] == helpers.OUTLIER_DO_NOT_USE
        result.shelve(zeroth, modify=False)

    # Verify that intermediate files are saved at the specified location
    i2d_files = glob(os.path.join(tmp_cwd, "*i2d.fits"))
    median_files = glob(os.path.join(tmp_cwd, "*median.fits"))
    assert len(i2d_files) != 0
    assert len(median_files) != 0


def test_outlier_step_on_disk(three_sci_as_asn, tmp_cwd):
    """Test whole step with an outlier including saving intermediate and results files"""
    container = ModelLibrary(three_sci_as_asn, on_disk=True)

    # Save all the data into a separate array before passing into step
    data_as_cube = list(
        container.map_function(lambda model, index: model.data.copy(), modify=False)
    )

    result = OutlierDetectionStep.call(
        container, save_results=True, save_intermediate_results=True, in_memory=False
    )

    with result:
        for i, r in enumerate(result):
            # Make sure nothing changed in SCI array except outliers are NaN
            dnu = (r.dq & helpers.OUTLIER_DO_NOT_USE).astype(bool)
            assert np.all(np.isnan(r.data[dnu]))
            assert_allclose(data_as_cube[i][~dnu], r.data[~dnu])

            # Verify source is not flagged
            assert r.dq[7, 7] == datamodels.dqflags.pixel["GOOD"]
            result.shelve(r, modify=False)

    # Verify CR is flagged
    with result:
        zeroth = result.borrow(0)
        assert zeroth.dq[12, 12] == helpers.OUTLIER_DO_NOT_USE
        result.shelve(zeroth, modify=False)

    # Verify intermediate results were written to disk
    dirname = tmp_cwd
    all_files = glob(os.path.join(dirname, "*.fits"))
    input_files = glob(os.path.join(dirname, "*_cal.fits"))
    result_files = glob(os.path.join(dirname, "*outlierdetectionstep.fits"))
    i2d_files = glob(os.path.join(dirname, "*i2d*.fits"))
    s2d_files = glob(os.path.join(dirname, "*outlier_s2d.fits"))
    median_files = glob(os.path.join(dirname, "*median.fits"))
    blot_files = glob(os.path.join(dirname, "*blot.fits"))

    assert len(result_files) == len(container)

    # i2d, median, blot files are written to the output directory
    assert len(i2d_files) == len(container)
    assert len(blot_files) == len(container)
    assert len(median_files) == 1

    # s2d files not written
    assert len(s2d_files) == 0

    # nothing else was written
    assert len(all_files) == len(input_files) + len(i2d_files) + len(median_files) + len(
        result_files
    ) + len(blot_files)


@pytest.mark.xfail(reason="Test data needs to be fixed to avoid outliers being detected.")
def test_outlier_step_square_source_no_outliers(mirimage_three_sci, tmp_cwd):
    """Test whole step with square source with sharp edges, no outliers"""
    container = ModelLibrary(list(mirimage_three_sci))

    # put a square source in all three exposures
    with container:
        for ccont in container:
            ccont.data[5:15, 5:15] += 1e3
            container.shelve(ccont)

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    dq_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data.copy())
            dq_as_cube.append(model.dq.copy())
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(container, in_memory=True)

    # Make sure nothing changed in SCI and DQ arrays
    with container:
        for i, image in enumerate(container):
            np.testing.assert_allclose(image.data, data_as_cube[i])
            np.testing.assert_allclose(image.dq, dq_as_cube[i])
            container.shelve(image, modify=False)

    # Make sure nothing changed in SCI and DQ arrays
    with result:
        for i, corrected in enumerate(result):
            np.testing.assert_allclose(data_as_cube[i], corrected.data)
            np.testing.assert_allclose(dq_as_cube[i], corrected.dq)
            result.shelve(corrected, modify=False)


def test_skip_one_exposure_imaging():
    model = datamodels.ImageModel()
    input_models = [model]
    step = OutlierDetectionStep()
    step.mode = "imaging"

    # Run the specified mode on one input model:
    # it should skip and return a copy of the input with the status set
    result = step.run(input_models)

    # Step is skipped
    assert query_step_status(result, "outlier_detection") == "SKIPPED"

    # Input is not modified
    assert result is not model
    assert model.meta.cal_step.outlier_detection is None
