import pytest
import numpy as np
from scipy.ndimage.filters import gaussian_filter

from jwst.outlier_detection.outlier_detection import flag_cr
from jwst import datamodels


OUTLIER_DO_NOT_USE = np.bitwise_or(datamodels.dqflags.pixel["DO_NOT_USE"],
                                    datamodels.dqflags.pixel["OUTLIER"])


@pytest.fixture
def sci_blot_image_pair():
    """Provide a science and blotted ImageModel pair."""
    shape = (10, 10)
    sigma = 0.02
    background = 3

    sci = datamodels.ImageModel(shape)

    # Populate keywords
    sci.meta.exposure.exposure_time = 1

    sci.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci.err = np.zeros(shape) + sigma

    # Add a source in the center
    sci.data[5, 5] += 20 * sigma

    # The blot image is just a smoothed version of the science image
    blot = sci.copy()
    blot.data = gaussian_filter(blot.data, sigma=3)

    return sci, blot


def test_flag_cr(sci_blot_image_pair):
    """Test the flag_cr function.  Test logic, not the actual noise model."""
    sci, blot = sci_blot_image_pair
    assert (sci.dq == 0).all()

    # Drop some CRs on the science array
    sci.data[3, 3] += 10
    sci.data[3, 7] += 100
    sci.data[7, 3] += 1e3
    sci.data[7, 7] += 1e4

    # run flag_cr() which updates in-place.  Copy sci first.
    data_copy = sci.data.copy()
    flag_cr(sci, blot)

    # Make sure science data array is unchanged after flag_cr()
    np.testing.assert_allclose(sci.data, data_copy)

    # Verify that both DQ flags are set in the DQ array for all outliers
    assert sci.dq[3, 3] == OUTLIER_DO_NOT_USE
    assert sci.dq[3, 7] == OUTLIER_DO_NOT_USE
    assert sci.dq[7, 3] == OUTLIER_DO_NOT_USE
    assert sci.dq[7, 7] == OUTLIER_DO_NOT_USE

    # Verify the source wasn't flagged
    assert sci.dq[5, 5] == datamodels.dqflags.pixel["GOOD"]


def test_flag_cr_with_subtracted_background(sci_blot_image_pair):
    """Test the flag_cr function on background-subtracted data"""
    sci, blot = sci_blot_image_pair

    sci.meta.background.subtracted = True
    sci.meta.background.level = 3

    # Subtract off the backgrounds
    sci.data -= 3
    blot.data -= 3

    # Drop a CR on the science array
    sci.data[8, 8] += 10

    # run flag_cr() which updates in-place.  Copy sci first.
    data_copy = sci.data.copy()
    flag_cr(sci, blot)

    # Make sure science data array is unchanged after flag_cr()
    np.testing.assert_allclose(sci.data, data_copy)

    # Verify that both DQ flags are set in the DQ array for all outliers
    assert sci.dq[8, 8] == OUTLIER_DO_NOT_USE

    # Verify the source wasn't flagged
    assert sci.dq[5, 5] == datamodels.dqflags.pixel["GOOD"]
