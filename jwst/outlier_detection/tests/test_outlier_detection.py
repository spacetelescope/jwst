import pytest
import numpy as np
from scipy.ndimage.filters import gaussian_filter

from jwst.outlier_detection import OutlierDetectionStep
from jwst.outlier_detection.outlier_detection import flag_cr
from jwst import datamodels
from jwst.assign_wcs.pointing import create_fitswcs


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


@pytest.fixture
def we_three_sci():
    """Provide 3 science images with different noise but identical source"""
    shape = (10, 10)
    sigma = 0.02
    background = 3

    sci = datamodels.ImageModel(shape)

    # Populate keywords
    sci.meta.instrument.name = "MIRI"
    sci.meta.instrument.detector = "MIRIMAGE"
    sci.meta.exposure.type = "MIR_IMAGE"
    sci.meta.observation.date = "2020-01-01"
    sci.meta.observation.time = "00:00:00"
    sci.meta.telescope = "JWST"
    sci.meta.exposure.exposure_time = 1
    sci.meta.wcsinfo.wcsaxes = 2
    sci.meta.wcsinfo.ctype1 = "RA---TAN"
    sci.meta.wcsinfo.ctype2 = "DEC--TAN"
    sci.meta.wcsinfo.cdelt1 = 1
    sci.meta.wcsinfo.cdelt2 = 1
    sci.meta.wcsinfo.roll_ref = 0
    sci.meta.wcsinfo.v3yangle = 0
    sci.meta.wcsinfo.vparity = -1
    sci.meta.wcsinfo.pc1_1 = -1
    sci.meta.wcsinfo.pc1_2 = 0
    sci.meta.wcsinfo.pc2_1 = 0
    sci.meta.wcsinfo.pc2_2 = 1
    sci.meta.wcsinfo.crpix1 = 5
    sci.meta.wcsinfo.crpix2 = 5
    sci.meta.wcsinfo.cunit1 = "deg"
    sci.meta.wcsinfo.cunit2 = "deg"
    sci.meta.filename = "foo.fits"

    sci.meta.wcs = create_fitswcs(sci)

    sci.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci.err = np.zeros(shape) + sigma

    # Add a source in the center
    sci.data[5, 5] += 20 * sigma

    # Make copies with different noise
    sci2 = sci.copy()
    sci3 = sci.copy()
    sci2.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci3.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci2.data[5, 5] += 20 * sigma
    sci3.data[5, 5] += 20 * sigma

    return sci, sci2, sci3


def test_outlier_step(we_three_sci):
    """Test whole step"""
    container = datamodels.ModelContainer(list(we_three_sci))
    result = OutlierDetectionStep.call(container)

    # Make sure nothing changed in SCI and DQ arrays
    for image, corrected in zip(container, result):
        np.testing.assert_allclose(image.data, corrected.data)
        np.testing.assert_allclose(image.dq, corrected.dq)

    # Drop some CRs on the science array
    container[0].data[3, 7] += 1e5
    container[1].data[7, 3] += 1e4
    container[1].data[7, 7] += 100
    container[2].data[3, 3] += 15

    # Run outlier again
    result = OutlierDetectionStep.call(container)

    # Verify CRs are flagged
    assert result[0].dq[3, 7] == OUTLIER_DO_NOT_USE
    assert result[1].dq[7, 3] == OUTLIER_DO_NOT_USE
    assert result[1].dq[7, 7] == OUTLIER_DO_NOT_USE
    assert result[2].dq[3, 3] == OUTLIER_DO_NOT_USE

    # Verify the SCI data remains unchanged
    for images, corrected in zip(container, result):
        np.testing.assert_allclose(images.data, corrected.data)
