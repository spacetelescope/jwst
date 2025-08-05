import os
from glob import glob

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ModelContainer
from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.tests import helpers
from jwst.resample.tests.test_resample_step import miri_rate_model


@pytest.mark.parametrize("resample", [True, False])
@pytest.mark.parametrize("save_intermediate", [True, False])
def test_outlier_step_spec(tmp_cwd, tmp_path, resample, save_intermediate):
    """Test outlier step for spec data including saving intermediate results."""
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)
    output_dir = str(output_dir)

    # Make a MIRI model and assign a spectral wcs
    miri_rate = miri_rate_model()
    miri_cal = AssignWcsStep.call(miri_rate)

    # Make it an exposure type outlier detection expects
    miri_cal.meta.exposure.type = "MIR_LRS-FIXEDSLIT"

    # Make a couple copies, give them unique exposure numbers and filename
    container = ModelContainer([miri_cal.copy(), miri_cal.copy(), miri_cal.copy()])
    for i, model in enumerate(container):
        model.meta.filename = f"test_{i}_cal.fits"

    # Drop a CR on the science array in the first image
    container[0].data[209, 37] += 1

    # Call outlier detection
    result = OutlierDetectionStep.call(
        container,
        resample_data=resample,
        output_dir=output_dir,
        save_results=True,
        save_intermediate_results=save_intermediate,
    )

    # Make sure nothing changed in SCI array
    for i, image in enumerate(result):
        nn = ~np.isnan(image.data)
        np.testing.assert_allclose(image.data[nn], miri_cal.data[nn])

        # Step is complete; input is not modified
        assert image.meta.cal_step.outlier_detection == "COMPLETE"
        assert image is not container[i]
        assert container[i].meta.cal_step.outlier_detection is None

    # Verify CR is flagged
    assert np.isnan(result[0].data[209, 37])
    assert result[0].dq[209, 37] == helpers.OUTLIER_DO_NOT_USE

    # Verify that intermediate files are saved at the specified location
    if save_intermediate:
        expected_intermediate = len(container)
    else:
        expected_intermediate = 0
    for dirname in [output_dir, tmp_cwd]:
        all_files = glob(os.path.join(dirname, "*.fits"))
        result_files = glob(os.path.join(dirname, "*outlierdetectionstep.fits"))
        i2d_files = glob(os.path.join(dirname, "*i2d*.fits"))
        s2d_files = glob(os.path.join(dirname, "*outlier_s2d.fits"))
        median_files = glob(os.path.join(dirname, "*median.fits"))
        blot_files = glob(os.path.join(dirname, "*blot.fits"))
        if dirname == output_dir:
            # Result files are always written to the output directory
            assert len(result_files) == len(container)

            # s2d and blot files are written to the output directory
            # if save_intermediate is True and resampling is set
            if resample:
                assert len(s2d_files) == expected_intermediate
                assert len(blot_files) == expected_intermediate
            else:
                assert len(s2d_files) == 0
                assert len(blot_files) == 0

            # Only one median file is saved if save_intermediate is True,
            # no matter how many input files there are
            if save_intermediate:
                assert len(median_files) == 1
            else:
                assert len(median_files) == 0

            # i2d files are never written
            assert len(i2d_files) == 0

            # Nothing else was written
            assert len(all_files) == (
                len(s2d_files) + len(median_files) + len(result_files) + len(blot_files)
            )
        else:
            # Nothing should be written to the current directory
            assert len(result_files) == 0
            assert len(s2d_files) == 0
            assert len(median_files) == 0
            assert len(i2d_files) == 0
            assert len(blot_files) == 0
            assert len(all_files) == 0

    miri_rate.close()
    result.close()
    for model in container:
        model.close()


def test_skip_one_exposure_spec():
    model = datamodels.ImageModel()
    input_models = [model]
    step = OutlierDetectionStep()
    step.mode = "spec"

    # Run the specified mode on one input model:
    # it should skip and return a copy of the input with the status set
    result = step.run(input_models)

    # Step is skipped
    assert result[0].meta.cal_step.outlier_detection == "SKIPPED"

    # Input is not modified
    assert result[0] is not model
    assert model.meta.cal_step.outlier_detection is None
