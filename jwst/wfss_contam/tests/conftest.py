import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.convolution import convolve
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils.datasets import make_100gaussians_image
from photutils.segmentation import SourceFinder, make_2dgaussian_kernel

from jwst.assign_wcs.tests.test_niriss import create_imaging_wcs, create_wfss_wcs

DIR_IMAGE = "direct_image.fits"


@pytest.fixture(scope="module")
def direct_image():
    """
    Make a mock direct image for testing.

    Returns
    -------
    np.ndarray
        Mock 2-d direct image data.
    """
    data = make_100gaussians_image(noise=False)
    kernel = make_2dgaussian_kernel(3, size=5)
    data = convolve(data, kernel)
    return data


@pytest.fixture(scope="module")
def direct_image_with_gradient(tmp_cwd_module, direct_image):  # noqa: ARG001
    """
    Add a gradient to the direct image and save it as a JWST datamodel.

    Returns
    -------
    ImageModel
        Direct image with a background gradient.
    """
    ny, nx = direct_image.shape
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / 5000.0
    data = direct_image + gradient

    # obs expects input list of direct image filenames
    model = dm.ImageModel(data=data)
    model.save(DIR_IMAGE)

    return model


@pytest.fixture(scope="module")
def segmentation_map(direct_image):
    """
    Make a segmentation map from the mock direct image.

    Returns
    -------
    SegmentationMapModel
        The segmentation map as a JwstDataModel.
    """
    _mean, median, stddev = sigma_clipped_stats(direct_image, sigma=3.0)
    threshold = median + 3 * stddev
    finder = SourceFinder(npixels=10)
    segm = finder(direct_image, threshold)

    # turn this into a jwst datamodel
    model = dm.SegmentationMapModel(data=segm.data)
    model.meta.wcs = create_imaging_wcs("F200W")
    return model


@pytest.fixture(scope="module")
def source_catalog(segmentation_map):
    """
    Create a mock source catalog from the segmentation map.

    Returns
    -------
    astropy.table.Table
        The source catalog with mock data.
    """
    source_ids = np.unique(segmentation_map.data)
    source_ids = source_ids[source_ids > 0]  # suppress background

    rng = np.random.default_rng(42)
    data = {
        "source_id": source_ids,
        "xcentroid": rng.uniform(0, segmentation_map.data.shape[1], len(source_ids)),
        "ycentroid": rng.uniform(0, segmentation_map.data.shape[0], len(source_ids)),
        "flux": rng.uniform(100, 1000, len(source_ids)),
    }
    return Table(data)


@pytest.fixture(scope="module")
def grism_wcs():
    """
    Load a WCS object from a grism model.

    Returns
    -------
    gwcs.wcs.WCS
        The grism wcs object.
    """
    return create_wfss_wcs("GR150C", pupil="F200W")
