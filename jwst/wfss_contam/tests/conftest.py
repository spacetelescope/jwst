import types

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.convolution import convolve
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils.datasets import make_100gaussians_image
from photutils.segmentation import SourceFinder, make_2dgaussian_kernel

from jwst.assign_wcs.tests.test_niriss import create_imaging_wcs, create_wfss_wcs


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
def direct_image_with_gradient(direct_image):
    """
    Add a gradient to the direct image and save it as a JWST datamodel.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.ImageModel`
        Direct image with a background gradient.
    """
    ny, nx = direct_image.shape
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / 5000.0
    data = direct_image + gradient

    # obs expects input list of direct image filenames
    model = dm.ImageModel(data=data)
    model.meta.wcs = create_imaging_wcs("F200W")

    return model


@pytest.fixture(scope="module")
def direct_image_cube_with_gradient(direct_image):
    """
    Build a multi-band direct image cube and save it as a WFSSCubeModel.

    Each wavelength plane is a copy of the direct image with a different linear
    gradient added, so the planes are not identical and tests can distinguish
    which wavelength was sampled.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.WFSSCubeModel`
        Multi-band direct image, also saved to filename "direct_image_cube.fits"
    """
    ny, nx = direct_image.shape
    y, x = np.mgrid[:ny, :nx]

    band_wls = np.array([1.75, 1.90, 2.05, 2.20], dtype=np.float32)
    n_bands = len(band_wls)
    cube = np.zeros((n_bands, ny, nx), dtype=np.float32)
    for i, wl in enumerate(band_wls):
        # Scale the gradient by wavelength so each plane is distinct
        gradient = wl * x * y / 5000.0
        cube[i] = direct_image + gradient

    model = dm.WFSSCubeModel(data=cube, wavelength=band_wls)
    model.meta.wcs = create_imaging_wcs("F200W")

    return model


@pytest.fixture(scope="module")
def segmentation_map(direct_image):
    """
    Make a segmentation map from the mock direct image.

    Yields
    ------
    SegmentationMapModel
        The segmentation map as a JwstDataModel.
    """
    _mean, median, stddev = sigma_clipped_stats(direct_image, sigma=3.0)
    threshold = median + 3 * stddev
    finder = SourceFinder(n_pixels=10)
    segm = finder(direct_image, threshold)

    # turn this into a jwst datamodel
    model = dm.SegmentationMapModel(data=segm.data)
    model.meta.wcs = create_imaging_wcs("F200W")
    yield model
    model.close()


def _sky_bbox(xcentroid, ycentroid, wcs, half_size=10):
    pix_bbox_ll = np.column_stack([xcentroid - half_size, ycentroid - half_size])
    pix_bbox_lr = np.column_stack([xcentroid + half_size, ycentroid - half_size])
    pix_bbox_ur = np.column_stack([xcentroid + half_size, ycentroid + half_size])
    pix_bbox_ul = np.column_stack([xcentroid - half_size, ycentroid + half_size])
    sky_bbox_ll = wcs.pixel_to_world(pix_bbox_ll[:, 0], pix_bbox_ll[:, 1])
    sky_bbox_lr = wcs.pixel_to_world(pix_bbox_lr[:, 0], pix_bbox_lr[:, 1])
    sky_bbox_ur = wcs.pixel_to_world(pix_bbox_ur[:, 0], pix_bbox_ur[:, 1])
    sky_bbox_ul = wcs.pixel_to_world(pix_bbox_ul[:, 0], pix_bbox_ul[:, 1])

    return sky_bbox_ll, sky_bbox_lr, sky_bbox_ur, sky_bbox_ul


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
    xcentroid = rng.uniform(0, segmentation_map.data.shape[1], len(source_ids))
    ycentroid = rng.uniform(0, segmentation_map.data.shape[0], len(source_ids))
    sky_bbox = _sky_bbox(xcentroid, ycentroid, segmentation_map.meta.wcs)
    data = {
        "label": source_ids,
        "xcentroid": xcentroid,
        "ycentroid": ycentroid,
        "isophotal_abmag": rng.uniform(20, 30, len(source_ids)),
        "sky_bbox_ll": sky_bbox[0],
        "sky_bbox_lr": sky_bbox[1],
        "sky_bbox_ur": sky_bbox[2],
        "sky_bbox_ul": sky_bbox[3],
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


@pytest.fixture(scope="module")
def phot_table():
    """
    Make a mock photom reference model for NIRISS WFSS.

    Returns
    -------
    :class:`~NisWfssPhotomModel`
        The photom reference model.
    """
    filt = [
        "GR150C",
    ] * 5
    pupil = [
        "F200W",
    ] * 5
    order = [-1, 0, 1, 2, 3]
    photmjsr = [1230.1855, 195.88435, 49.010494, 1584.0671, 7298.0225]
    uncertainty = [12.301856, 1.9588436, 0.49010494, 15.840671, 72.980225]
    nelem = [100] * 5
    wavelength = [
        np.linspace(1.7, 2.3, 100),
    ] * 5

    response_scaling = [10, 1, 1, 10, 100]  # response is inversely proportional to throughput
    relresponse = [np.ones_like(wavelength[0]) * r for r in response_scaling]
    reluncertainty = [np.ones_like(wavelength[0]) * r * 0.01 for r in response_scaling]

    dtype = [
        ("filter", "S12"),
        ("pupil", "S15"),
        ("order", "<i2"),
        ("photmjsr", "<f4"),
        ("uncertainty", "<f4"),
        ("nelem", "<i2"),
        ("wavelength", "<f4", (wavelength[0].size,)),
        ("relresponse", "<f4", (relresponse[0].size,)),
        ("reluncertainty", "<f4", (reluncertainty[0].size,)),
    ]
    phot = np.recarray((5,), dtype=dtype)

    phot["filter"] = filt
    phot["pupil"] = pupil
    phot["order"] = order
    phot["photmjsr"] = photmjsr
    phot["uncertainty"] = uncertainty
    phot["nelem"] = nelem
    phot["wavelength"] = wavelength
    phot["relresponse"] = relresponse
    phot["reluncertainty"] = reluncertainty
    return phot


@pytest.fixture(scope="module")
def photom_ref_model_niriss(phot_table):
    """
    Make a mock photom reference model for NIRISS WFSS.

    Parameters
    ----------
    phot_table : np.recarray
        Photometry table.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.NisWfssPhotomModel`
        Photom ref file model.
    """
    model = dm.NisWfssPhotomModel(phot_table=phot_table)
    model.phot_unit = "MJy micron s / (DN sr)"
    yield model
    model.close()


@pytest.fixture(scope="module")
def photom_ref_model_nircam(phot_table):
    """
    Make a mock photom reference model for NIRCam WFSS.

    Parameters
    ----------
    phot_table : np.recarray
        Photometry table.

    Yields
    ------
    `~stdatamodels.jwst.datamodels.NrcWfssPhotomModel`
        Photom ref file model.
    """
    model = dm.NrcWfssPhotomModel(phot_table=phot_table)
    model.phot_unit = "MJy Angstrom s / (DN sr)"
    yield model
    model.close()


@pytest.fixture
def photom_ref_model(request):
    """
    Make a mock photom reference model that can be parameterized for any WFSS mode.

    If request.param is not set, defaults to NIRISS.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.JwstDataModel`
        Photom ref file model of the appropriate type based on the parameterization.
    """
    fixture_name = getattr(request, "param", "photom_ref_model_niriss")
    return request.getfixturevalue(fixture_name)


@pytest.fixture(scope="module")
def wavelengthrange_ref_model():
    """
    Mock WAVELENGTHRANGE reference file object simulating NIRISS order 1.

    Returns
    -------
    types.SimpleNamespace
        Mock class with `order` and `get_wfss_wavelength_range` attributes.
    """
    wr = types.SimpleNamespace()
    wr.order = np.array([1])
    wr.get_wfss_wavelength_range = lambda _filter_name, orders: {orders[0]: (1.0, 3.0)}
    return wr
