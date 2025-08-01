import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from stdatamodels.jwst import datamodels

from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.tests import helpers


def test_outlier_step_weak_cr_coron(miri_container, tmp_cwd):
    """Test whole step with an outlier for an example coronagraphic mode"""
    container = miri_container

    # Drop a weak CR on the science array
    # no noise so it should always be above the default threshold of 5
    container[0].data[12, 12] = helpers.BACKGROUND + helpers.SIGMA * 10

    # coron3 will provide a CubeModel so convert the container to a cube
    cube = helpers.container_to_cube(container)

    result = OutlierDetectionStep.call(cube)

    # Step is complete; input is not modified
    assert result.meta.cal_step.outlier_detection == "COMPLETE"
    assert result is not cube
    assert cube.meta.cal_step.outlier_detection is None

    # Make sure nothing changed in SCI array except that
    # outliers are NaN
    for i, image in enumerate(container):
        dnu = (result.dq[i] & helpers.OUTLIER_DO_NOT_USE).astype(bool)
        assert np.all(np.isnan(result.data[i][dnu]))
        assert_allclose(image.data[~dnu], result.data[i][~dnu])

    # Verify source is not flagged
    j, k = helpers.SIGNAL_LOC
    assert_array_equal(result.dq[:, j, k], datamodels.dqflags.pixel["GOOD"])

    # Verify CR is flagged
    assert result.dq[0, 12, 12] == helpers.OUTLIER_DO_NOT_USE
    assert np.isnan(result.data[0, 12, 12])


def test_coron_save_intermediate(miri_container, tmp_path):
    """Test intermediate files for coronagraphic mode"""
    container = miri_container
    cube = helpers.container_to_cube(container)

    OutlierDetectionStep.call(
        cube, output_dir=str(tmp_path), save_results=True, save_intermediate_results=True
    )
    basename = "foo1"
    assert (tmp_path / f"{basename}_median.fits").exists()
    assert (tmp_path / f"{basename}_outlierdetectionstep.fits").exists()


def test_coron_error_invalid():
    model = datamodels.ImageModel()
    step = OutlierDetectionStep()
    step.mode = "coron"

    # model containers not allowed
    with pytest.raises(TypeError, match="must be a CubeModel"):
        step.run(model)
