import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.tests import helpers


@pytest.mark.parametrize("rolling_window_width", [7, 0])
def test_outlier_step_weak_cr_tso(mirimage_50_sci, rolling_window_width):
    """Test outlier detection with rolling median on time-varying source
    This test fails if rolling_window_width is set to 0, i.e., take simple median
    """
    im = ModelContainer(mirimage_50_sci)

    # Drop a weak CR on the science array
    cr_timestep = 5
    im[cr_timestep].data[12, 12] = helpers.BACKGROUND + helpers.SIGMA * 10

    # make time variability that has larger total amplitude than
    # the CR signal but deviations frame-by-frame are smaller
    j, k = helpers.SIGNAL_LOC
    real_time_variability = helpers.SIGNAL * np.cos(np.linspace(0, np.pi, 50))
    for i, model in enumerate(im):
        model.data[j, k] += real_time_variability[i]
        model.err[j, k] = np.sqrt(helpers.SIGMA**2 + model.data[j, k])

    cube = helpers.container_to_cube(im)

    result = OutlierDetectionStep.call(cube, rolling_window_width=rolling_window_width)

    # Step is complete; input is not modified
    assert result.meta.cal_step.outlier_detection == "COMPLETE"
    assert result is not cube
    assert cube.meta.cal_step.outlier_detection is None

    # Make sure nothing changed in SCI array except
    # that outliers are NaN
    for i, model in enumerate(im):
        dnu = (result.dq[i] & helpers.OUTLIER_DO_NOT_USE).astype(bool)
        assert np.all(np.isnan(result.data[i][dnu]))
        assert_allclose(model.data[~dnu], result.data[i][~dnu])

    # Verify source is not flagged for rolling median
    if rolling_window_width == 7:
        assert_array_equal(result.dq[:, j, k], datamodels.dqflags.pixel["GOOD"])
    # But this fails for simple median
    elif rolling_window_width == 0:
        with pytest.raises(AssertionError):
            assert_array_equal(result.dq[:, j, k], datamodels.dqflags.pixel["GOOD"])

    # Verify CR is flagged
    assert result.dq[cr_timestep, 12, 12] == helpers.OUTLIER_DO_NOT_USE


def test_tso_save_intermediate(tmp_path, mirimage_50_sci):
    """Test intermediate output for TSO."""
    im = ModelContainer(mirimage_50_sci)
    cube = helpers.container_to_cube(im)

    OutlierDetectionStep.call(
        cube,
        output_dir=str(tmp_path),
        save_results=True,
        save_intermediate_results=True,
    )
    basename = "foo1"
    assert (tmp_path / f"{basename}_median.fits").exists()
    assert (tmp_path / f"{basename}_outlierdetectionstep.fits").exists()


def test_tso_error_invalid():
    model = ModelContainer([datamodels.CubeModel()])
    step = OutlierDetectionStep()
    step.mode = "tso"

    # model containers not allowed
    with pytest.raises(TypeError, match="does not support ModelContainer input"):
        step.run(model)
