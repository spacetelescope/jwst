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
    shape = (20, 20)
    sigma = 0.02
    background = 3

    sci = datamodels.ImageModel(shape)

    # Populate keywords
    sci.meta.exposure.exposure_time = 1
    sci.meta.background.subtracted = False
    sci.meta.background.level = background

    sci.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci.err = np.zeros(shape) + sigma

    # Add a source in the center
    sci.data[10, 10] += 20 * sigma

    # The blot image is just a smoothed version of the science image that has
    # its background subtracted
    blot = sci.copy()
    blot.data = gaussian_filter(blot.data, sigma=3)
    blot.data -= background

    return sci, blot


def test_flag_cr(sci_blot_image_pair):
    """Test the flag_cr function.  Test logic, not the actual noise model."""
    sci, blot = sci_blot_image_pair
    assert (sci.dq == 0).all()

    # Drop some CRs on the science array
    sci.data[3, 3] += 100
    sci.data[3, 7] += 1e3
    sci.data[7, 3] += 1e4
    sci.data[7, 7] += 1e5

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
    assert sci.dq[10, 10] == datamodels.dqflags.pixel["GOOD"]


@pytest.fixture
def we_three_sci():
    """Provide 3 science images with different noise but identical source
    and same background level"""
    shape = (20, 20)
    sigma = 0.02
    background = 1.5
    signal = 200 * sigma

    sci1 = datamodels.ImageModel(shape)

    # Populate keywords
    sci1.meta.instrument.name = "MIRI"
    sci1.meta.instrument.detector = "MIRIMAGE"
    sci1.meta.exposure.type = "MIR_IMAGE"
    sci1.meta.observation.date = "2020-01-01"
    sci1.meta.observation.time = "00:00:00"
    sci1.meta.telescope = "JWST"
    sci1.meta.exposure.exposure_time = 1
    sci1.meta.wcsinfo.wcsaxes = 2
    sci1.meta.wcsinfo.ctype1 = "RA---TAN"
    sci1.meta.wcsinfo.ctype2 = "DEC--TAN"
    sci1.meta.wcsinfo.cdelt1 = 3e-6
    sci1.meta.wcsinfo.cdelt2 = 3e-6
    sci1.meta.wcsinfo.roll_ref = 0
    sci1.meta.wcsinfo.v3yangle = 0
    sci1.meta.wcsinfo.vparity = -1
    sci1.meta.wcsinfo.pc1_1 = 1
    sci1.meta.wcsinfo.pc1_2 = 0
    sci1.meta.wcsinfo.pc2_1 = 0
    sci1.meta.wcsinfo.pc2_2 = 1
    sci1.meta.wcsinfo.crpix1 = 5
    sci1.meta.wcsinfo.crpix2 = 5
    sci1.meta.wcsinfo.crval1 = 0
    sci1.meta.wcsinfo.crval2 = 0
    sci1.meta.wcsinfo.cunit1 = "deg"
    sci1.meta.wcsinfo.cunit2 = "deg"
    sci1.meta.background.subtracted = False
    sci1.meta.background.level = background

    # Replace the FITS-type WCS with an Identity WCS
    sci1.meta.wcs = create_fitswcs(sci1)


    sci1.err = np.zeros(shape) + sigma

    # Make copies with different noise
    sci2 = sci1.copy()
    sci3 = sci1.copy()

    # Populate data background with random noise
    sci1.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci2.data = np.random.normal(loc=background, size=shape, scale=sigma)
    sci3.data = np.random.normal(loc=background, size=shape, scale=sigma)

    sci1.meta.filename = "foo1_cal.fits"
    sci2.meta.filename = "foo2_cal.fits"
    sci3.meta.filename = "foo3_cal.fits"

    # Add a source in the centers
    sci1.data[7, 7] += signal
    sci2.data[7, 7] += signal
    sci3.data[7, 7] += signal

    return sci1, sci2, sci3


def test_outlier_step_no_outliers(we_three_sci):
    """Test whole step, no outliers"""
    container = datamodels.ModelContainer(list(we_three_sci))
    pristine = container.copy()
    result = OutlierDetectionStep.call(container)

    # Make sure nothing changed in SCI and DQ arrays
    for image, uncorrected in zip(pristine, container):
        np.testing.assert_allclose(image.data, uncorrected.data)
        np.testing.assert_allclose(image.dq, uncorrected.dq)

    # Make sure nothing changed in SCI and DQ arrays
    for image, corrected in zip(container, result):
        np.testing.assert_allclose(image.data, corrected.data)
        np.testing.assert_allclose(image.dq, corrected.dq)


def test_outlier_step(we_three_sci):
    """Test whole step with an outlier"""
    container = datamodels.ModelContainer(list(we_three_sci))

    # Drop a CR on the science array
    container[0].data[12, 12] += 1e3

    result = OutlierDetectionStep.call(container)

    # Make sure nothing changed in SCI array
    for image, corrected in zip(container, result):
        np.testing.assert_allclose(image.data, corrected.data)

    # Verify source is not flagged
    for r in result:
        assert r.dq[7, 7] == datamodels.dqflags.pixel["GOOD"]

    # Verify CR is flagged
    assert result[0].dq[12, 12] == OUTLIER_DO_NOT_USE
