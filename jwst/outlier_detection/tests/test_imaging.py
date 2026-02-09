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


@pytest.mark.parametrize("weight", ["exptime", "ivm"])
@pytest.mark.parametrize("src_type", ["square", "gaussian"])
def test_outlier_step_with_source_no_outliers(mirimage_three_sci, tmp_cwd, src_type, weight):
    """Test whole step with no outliers and an artificial source:
    a uniform "square" with sharp edges or a "gaussian" source centered
    on the image, no outliers."""
    container = ModelLibrary(list(mirimage_three_sci))

    atol = 2.0 * np.finfo(np.float32).eps

    # Create artificial source
    if src_type == "square":
        src = np.full((11, 11), 50 * helpers.SIGMA, dtype=np.float32)
    else:  # gaussian
        y, x = np.indices((11, 11)) - 5
        fwhm = 1.5
        src = 50 * helpers.SIGMA * np.exp(-4.0 * np.log(2.0) * (x**2 + y**2) / fwhm**2)

    # put a Gaussian source in all three exposures
    with container:
        for ccont in container:
            ccont.data[5:16, 5:16] += src
            ccont.err[5:16, 5:16] = np.sqrt(ccont.err[5:16, 5:16] ** 2 + src)
            container.shelve(ccont)

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    dq_as_cube = []
    non_nan_mask_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data.copy())
            dq_as_cube.append(model.dq.copy())
            non_nan_mask_as_cube.append(np.isfinite(model.data))
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(container, in_memory=True, weight_type=weight)

    # Make sure nothing changed in SCI and DQ arrays
    with container:
        for i, image in enumerate(container):
            m = non_nan_mask_as_cube[i]
            assert np.all(np.logical_and(np.isfinite(image.data), m))
            np.testing.assert_allclose(image.data[m], data_as_cube[i][m], rtol=0.0, atol=atol)
            np.testing.assert_array_equal(image.dq[m], dq_as_cube[i][m])
            container.shelve(image, modify=False)

    # Make sure nothing changed in SCI and DQ arrays
    with result:
        for i, corrected in enumerate(result):
            m = non_nan_mask_as_cube[i]
            assert np.all(np.logical_and(np.isfinite(corrected.data), m))
            np.testing.assert_allclose(data_as_cube[i][m], corrected.data[m], rtol=0.0, atol=atol)
            np.testing.assert_array_equal(dq_as_cube[i][m], corrected.dq[m])
            result.shelve(corrected, modify=False)


@pytest.mark.parametrize("weight", ["exptime", "ivm"])
@pytest.mark.parametrize("src_type", ["gaussian", "square", "strip", "point"])
@pytest.mark.parametrize("idx", [0, 1, 2])
def test_outlier_step_with_outliers(mirimage_three_sci, tmp_cwd, src_type, weight, idx):
    """Test whole step with outlier(s) added besides a uniform square source."""
    container = ModelLibrary(list(mirimage_three_sci))

    atol = 2.0 * np.finfo(np.float32).eps

    # Create artificial source
    src_sl = np.s_[5:16, 5:16]
    src = np.full((11, 11), 100 * helpers.SIGMA, dtype=np.float32)

    # Create artificial CR
    if src_type == "strip":
        sl = np.s_[9:12, 5:16]
        cr = np.full((3, 11), 100 * helpers.SIGMA, dtype=np.float32)
    elif src_type == "square":
        sl = np.s_[5:16, 5:16]
        cr = np.full((11, 11), 100 * helpers.SIGMA, dtype=np.float32)
    elif src_type == "point":
        sl = ((1, 1), (10, 10))
        cr = 100 * helpers.SIGMA
    else:  # gaussian
        sl = np.s_[5:16, 5:16]
        y, x = np.indices((11, 11)) - 5
        fwhm = 1.5
        cr = 100 * helpers.SIGMA * np.exp(-4.0 * np.log(2.0) * (x**2 + y**2) / fwhm**2)

    cr_mask = cr > 5 * helpers.SIGMA

    # put a Gaussian source in all three exposures
    with container:
        for i, ccont in enumerate(container):
            ccont.data[src_sl] += src
            if i == idx:
                ccont.data[sl] += cr
                # ccont.err[5:16, 5:16] = np.sqrt(ccont.err[5:16, 5:16]**2 + src)
            container.shelve(ccont)

    # Save all the data into a separate array before passing into step
    data_as_cube = []
    dq_as_cube = []
    non_nan_mask_as_cube = []
    with container:
        for model in container:
            data_as_cube.append(model.data.copy())
            dq_as_cube.append(model.dq.copy())
            non_nan_mask_as_cube.append(np.isfinite(model.data))
            container.shelve(model, modify=False)

    result = OutlierDetectionStep.call(container, in_memory=True, weight_type=weight)

    with result:
        for i, corrected in enumerate(result):
            m = non_nan_mask_as_cube[i]
            if i == idx:
                # Make sure CR pixels are now NaN in SCI and flagged in DQ
                m[sl] = np.logical_not(cr_mask)
                m2 = np.isfinite(corrected.data)

                if src_type == "square":
                    pytest.xfail(
                        "Square CR fails with pixels in the interior of the square not flagged as outliers."
                    )

                assert np.all(m == m2)
                assert np.all(
                    (corrected.dq[~m] & helpers.OUTLIER_DO_NOT_USE) == helpers.OUTLIER_DO_NOT_USE
                )

            np.testing.assert_allclose(data_as_cube[i][m], corrected.data[m], rtol=0.0, atol=atol)
            np.testing.assert_array_equal(dq_as_cube[i][m], corrected.dq[m])
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
