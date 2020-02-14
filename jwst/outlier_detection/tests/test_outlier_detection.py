import pytest
import numpy as np
from scipy.ndimage.filters import gaussian_filter

from jwst.outlier_detection.outlier_detection import flag_cr
from jwst import datamodels


@pytest.fixture
def sci_blot_image_pair():
    """Provide a science and blotted ImageModel pair."""
    shape = (10, 10)
    sci = datamodels.ImageModel(shape)

    # Populate keywords
    sci.meta.exposure.exposure_time = 1

    # Add poisson noise to image data
    p = np.random.poisson(size=shape, lam=1e3)
    sci.data = p / p.mean() - 1

    # The blot image is just a smoothed version of the science image
    blot = sci.copy()
    blot.data = gaussian_filter(blot.data, sigma=3)

    return sci, blot


def test_flag_cr(sci_blot_image_pair):
    """Test the flag_cr function.  Test logic, not the actual noise model."""
    sci, blot = sci_blot_image_pair
    assert (sci.dq == 0).all()

    # Add some background
    sci.data += 3
    blot.data += 3

    # Drop a CR on the science array
    sci.data[5, 5] += 10

    flag_cr(sci, blot)
    assert sci.dq[5, 5] > 0


def test_flag_cr_with_subtracted_background(sci_blot_image_pair):
    """Test the flag_cr function on background-subtracted data"""
    sci, blot = sci_blot_image_pair

    sci.meta.background.subtracted = True
    sci.meta.background.level = 3

    # Drop a CR on the science array
    sci.data[5, 5] += 10

    flag_cr(sci, blot)
    assert sci.dq[5, 5] > 0
