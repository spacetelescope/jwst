import warnings

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_miri import create_hdul
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file
from jwst.datamodels import ModelContainer
from jwst.outlier_detection.tests import helpers


@pytest.fixture()
def sci_blot_image_pair():
    """
    Provide a science and blotted ImageModel pair.

    Returns
    -------
    sci : ImageModel
        Science model.
    blot : ImageModel
        Blotted model.
    """
    background = 3

    sci = datamodels.ImageModel(helpers.SHAPE)

    # Populate keywords
    sci.meta.exposure.exposure_time = 1
    sci.meta.background.subtracted = False
    sci.meta.background.level = background

    rng = np.random.default_rng(720)
    sci.data = rng.normal(loc=background, size=helpers.SHAPE, scale=helpers.SIGMA)
    sci.err = np.zeros(helpers.SHAPE) + helpers.SIGMA
    sci.var_rnoise += 0

    # Add a source in the center
    signal = 20 * helpers.SIGMA
    sci.data[10, 10] += signal
    # update the noise for this source to include the photon/measurement noise
    sci.err[10, 10] = np.sqrt(helpers.SIGMA**2 + signal)

    # The blot image is just a smoothed version of the science image that has
    # its background subtracted
    blot = sci.copy()
    blot.data = gaussian_filter(blot.data, sigma=3)
    blot.data -= background

    return sci, blot


@pytest.fixture()
def scimodel_base():
    """
    Make a base model for all instruments and modes.

    Instrument-specific keywords should be modified by individual tests
    before AssignWCS is called.

    Returns
    -------
    ImageModel
        A model with basic metadata.
    """
    sci1 = datamodels.ImageModel(helpers.SHAPE)

    # Populate keywords
    sci1.meta.instrument.name = "MIRI"
    sci1.meta.instrument.detector = "MIRIMAGE"
    sci1.meta.exposure.type = "MIR_IMAGE"
    sci1.meta.visit.tsovisit = False
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
    sci1.meta.wcsinfo.ra_ref = 1.5e-5
    sci1.meta.wcsinfo.dec_ref = 1.5e-5
    sci1.meta.wcsinfo.v2_ref = 0
    sci1.meta.wcsinfo.v3_ref = 0
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
    sci1.meta.background.level = 1.5
    sci1.meta.target.ra = 0.0
    sci1.meta.target.dec = 0.0

    sci1.meta.filename = "foo1_cal.fits"

    # add pixel areas
    sci1.meta.photometry.pixelarea_steradians = 1.0
    sci1.meta.photometry.pixelarea_arcsecsq = 1.0

    # make some mock data with a "real" source at 7,7
    data, err = helpers.mock_data()
    sci1.data = data
    sci1.err = err
    sci1.var_rnoise = np.zeros(helpers.SHAPE) + 1.0

    return sci1


@pytest.fixture()
def three_sci(scimodel_base):
    """
    Make 3 science images.

    Each image will have different noise but identical source and background.

    Returns
    -------
    list of ImageModel
        3 science models.
    """
    all_sci = [scimodel_base]
    for i in range(2):
        tsci = scimodel_base.copy()
        data, err = helpers.mock_data()
        # Add a source in the center
        tsci.data = data
        tsci.err = err
        tsci.meta.filename = f"foo{i + 2}_cal.fits"
        all_sci.append(tsci)
    return all_sci


@pytest.fixture()
def mirimage_three_sci(three_sci):
    """
    Provide 3 MIRI imaging science observations.

    This is the default/base model set for the test suite, everything is default.
    So just need to assign the WCS, identical for all.
    This fixture is separated from three_sci so the latter
    can be reused for other instruments and modes

    Returns
    -------
    list of ImageModel
        3 science models with a WCS assigned.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=UserWarning, message="Double sampling check FAILED"
        )
        result = helpers.assign_wcs_to_models(three_sci, "MIR_IMAGE", False)
    return result


@pytest.fixture()
def mirimage_50_sci(scimodel_base):
    """
    Provide 50 MIRI TSO imaging science observations.

    Returns
    -------
    list of ImageModel
        50 science models.
    """
    # first call AssignWcsStep on the base model, all WCSs will be the same after copy
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=UserWarning, message="Double sampling check FAILED"
        )
        scimodel_base = AssignWcsStep.call(scimodel_base)
    scimodel_base.meta.visit.tsovisit = True

    all_sci = [scimodel_base]
    for i in range(49):
        tsci = scimodel_base.copy()
        data, err = helpers.mock_data()
        # Add a source in the center
        tsci.data = data
        tsci.err = err
        tsci.meta.filename = f"foo{i + 2}_cal.fits"
        all_sci.append(tsci)

    return all_sci


@pytest.fixture()
def three_sci_as_asn(mirimage_three_sci, tmp_cwd):  # noqa: ARG001
    """
    Create an association with 3 science images.

    The science images are saved to a temporary directory so
    the association can locate them when loaded.

    Returns
    -------
    dict
        ASN dictionary.
    """
    for model in mirimage_three_sci:
        model.save(model.meta.filename)
    filenames = [model.meta.filename for model in mirimage_three_sci]
    # Drop a CR on the science array
    with datamodels.open(filenames[0]) as dm0:
        dm0.data[12, 12] += 1
        dm0.save(dm0.meta.filename)

    # Initialize inputs for the test based on filenames only
    # needs to be an asn for ModelLibrary to load it in on_disk mode
    asn = {
        "asn_type": "test",
        "asn_id": "o001",
        "products": [
            {
                "name": "product_a",
                "members": [
                    {"expname": filenames[0], "exptype": "science"},
                    {"expname": filenames[1], "exptype": "science"},
                    {"expname": filenames[2], "exptype": "science"},
                ],
            },
        ],
    }
    return asn


@pytest.fixture()
def miri_container(three_sci):
    """
    Make a container from three science models, with a coron exp_type.

    Returns
    -------
    ModelContainer
        Container of 3 models with MIR_LYOT exp_type and an assigned WCS.
    """
    exptype = "MIR_LYOT"
    with pytest.warns(UserWarning, match="Double sampling check FAILED"):
        three_sci = helpers.assign_wcs_to_models(three_sci, exptype, False)
    container = ModelContainer(three_sci)
    return container


@pytest.fixture(scope="module")
def miri_ifu_rate():
    """
    Create a MIRI IFU test model.

    Test data is borrowed from assign_wcs tests.

    Yields
    ------
    IFUImageModel
        A MIRI IFU model with a wcs assigned.
    """
    # metadata for assign_wcs
    hdul = create_hdul(detector="MIRIFUSHORT", channel="12", band="SHORT")

    # mock data for processing
    shape = (100, 100)
    model = datamodels.IFUImageModel(hdul)
    model.data = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    model.err = np.full(shape, 1.0)
    model.dq = np.full(shape, 0)

    model_wcs = AssignWcsStep.call(model)
    hdul.close()
    model.close()
    yield model_wcs
    model_wcs.close()


@pytest.fixture(scope="module")
def nirspec_ifu_rate():
    """
    Create a NIRSpec IFU test model.

    Test data is borrowed from assign_wcs tests.

    Yields
    ------
    IFUImageModel
        A NIRSpec IFU model with a wcs assigned.
    """
    shape = (2048, 2048)
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    model = datamodels.IFUImageModel(hdul)
    model.data = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    model.err = np.full(shape, 1.0)
    model.dq = np.full(shape, 0)
    model_wcs = AssignWcsStep.call(model)

    hdul.close()
    model.close()
    yield model_wcs
    model_wcs.close()
