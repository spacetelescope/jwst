import json
from pathlib import Path

import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import models

# Constants for MIRI imager and LRS slit dimensions
PIXSCALE = 0.11  # arcsec/pixel for MIRI imager
JWST_DIAMETER = 6.5  # meters
X_REF = 326.13  # for these tests we assume same ref position for slit and slitless
Y_REF = 300.7
MIRI_DETECTOR_SHAPE = (1024, 1032)  # (ny, nx) for MIRI imager


def get_wavelength(filter_name):  # numpydoc ignore=RT01
    """Map filter name to central wavelength in microns."""
    filter_wavelengths = {
        "F560W": 5.6,
        "F770W": 7.7,
        "F1000W": 10.0,
        "F1130W": 11.3,
        "F1280W": 12.8,
        "F1500W": 15.0,
        "F1800W": 18.0,
        "F2100W": 21.0,
        "F2550W": 25.5,
    }
    return filter_wavelengths.get(filter_name, None)


def make_slitless_data(wavelength=15.0, offset=(0, 0)):
    """
    Make a fake MIRI LRS slitless TA verification image.

    Parameters
    ----------
    wavelength : float, optional
        Wavelength in microns.
    offset : tuple, optional
        (x,y) offset of the source from the reference center in pixels.

    Returns
    -------
    data : 2D ndarray
        Simulated slitless TA image.
    """
    # Create coordinate grids using full MIRI detector size
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]

    # Calculate diffraction-limited Airy disk radius
    # First zero of Airy disk: theta = 1.22 * lambda / D (in radians)
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slitless reference position
    x_center = X_REF + offset[0]
    y_center = Y_REF + offset[1]

    # Simulate source as an Airy disk model at diffraction-limited resolution
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data = airy(x, y)

    # Add some noise with fixed seed for reproducibility
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data.shape)
    data += noise

    return data


def make_ta_model(data):
    """
    Create a MIRI TA verification image model.

    Parameters
    ----------
    data : ndarray
        The image data array.

    Returns
    -------
    ImageModel
        MIRI TA verification image model with metadata.
    """
    model = dm.ImageModel()
    model.data = data
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "F1500W"
    model.meta.exposure.type = "MIR_TACONFIRM"

    model.meta.observation.date = "2025-11-10"
    model.meta.observation.time = "12:00:00"

    # Add subarray metadata
    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1  # 1-indexed FITS convention
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = MIRI_DETECTOR_SHAPE[1]  # nx
    model.meta.subarray.ysize = MIRI_DETECTOR_SHAPE[0]  # ny

    # add wcsinfo
    model.meta.wcsinfo.crpix1 = 0
    model.meta.wcsinfo.crpix2 = 0
    model.meta.wcsinfo.v2_ref = -414.848
    model.meta.wcsinfo.v3_ref = -400.426
    model.meta.wcsinfo.roll_ref = 92.69243893875218
    model.meta.wcsinfo.ra_ref = 108.63038309653392
    model.meta.wcsinfo.dec_ref = 13.860413558600543
    model.meta.wcsinfo.v3yangle = 4.75797
    model.meta.wcsinfo.vparity = -1
    model.meta.wcsinfo.crval1 = 108.63038309653392
    model.meta.wcsinfo.crval2 = 13.860413558600543

    # dither info
    model.meta.dither.x_offset = 0.0
    model.meta.dither.y_offset = 0.0

    return model


def make_slit_data(offset):
    """
    Make a fake MIRI LRS slit TA verification image.

    The output data are zero everywhere except in the slit region,
    where all data are weighted according to the pathloss correction.
    An Airy disk PSF is created at the source position plus offset,
    noise is added, and then the weights are applied.

    Parameters
    ----------
    offset : tuple
        (x,y) offset of the source from the reference center in pixels.

    Returns
    -------
    model : ImageModel
        Simulated slit TA image model with PSF weighted by pathloss.
    """
    wavelength = get_wavelength("F1500W")

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slit reference position
    x_center = X_REF + offset[0]
    y_center = Y_REF + offset[1]

    # Create Airy disk model on the full detector
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data_full = airy(x, y)

    # Make a slit mask
    slit_half_height = 8  # pixels in y direction (cross-dispersion) very approximate
    slit_half_width = 60  # pixels in x direction (along slit)
    y_min = max(0, int(Y_REF - slit_half_height))
    y_max = min(MIRI_DETECTOR_SHAPE[0], int(Y_REF + slit_half_height))
    x_min = max(0, int(X_REF - slit_half_width))
    x_max = min(MIRI_DETECTOR_SHAPE[1], int(X_REF + slit_half_width))
    slit_mask = np.zeros(MIRI_DETECTOR_SHAPE, dtype=bool)
    slit_mask[y_min:y_max, x_min:x_max] = True

    # Add noise
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data_full.shape)
    data_full += noise

    # Apply slit weights to the data
    data = data_full * slit_mask

    # Add bad values near the slit source position to test their handling
    data[int(Y_REF) + 4, int(X_REF) + 2] = np.nan
    data[int(Y_REF) - 3, int(X_REF) - 4] = np.inf

    model = make_ta_model(data)
    return model


def make_empty_lrs_model():
    """
    Create empty ImageModel with metadata shared between slit and slitless modes for LRS.

    Returns
    -------
    ImageModel
        MIRI LRS model with required metadata.
    """
    model = dm.ImageModel()
    model.meta.instrument.name = "MIRI"
    model.meta.target.source_type = "POINT"
    model.meta.instrument.filter = "F1500W"
    model.meta.instrument.detector = "MIRIMAGE"

    model.meta.observation.date = "2025-11-10"
    model.meta.observation.time = "12:00:00"

    # Add subarray metadata
    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1  # 1-indexed FITS convention
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = MIRI_DETECTOR_SHAPE[1]  # nx
    model.meta.subarray.ysize = MIRI_DETECTOR_SHAPE[0]  # ny

    # bare minimum wcs info to get assign_wcs step to pass
    model.meta.wcsinfo.crpix1 = 693.5
    model.meta.wcsinfo.crpix2 = 512.5
    model.meta.wcsinfo.v2_ref = -453.37849
    model.meta.wcsinfo.v3_ref = -373.810549
    model.meta.wcsinfo.roll_ref = 272.3237653262276
    model.meta.wcsinfo.ra_ref = 80.54724018120017
    model.meta.wcsinfo.dec_ref = -69.5081101864959
    model.meta.wcsinfo.v3yangle = 0.0
    model.meta.wcsinfo.vparity = 1

    return model


def make_ta_association(sci_model, ta_model=None, asn_fname="mir_lrs_ta_asn.json"):
    """
    Create an association file for TA centering test.

    Parameters
    ----------
    sci_model : ImageModel
        Science exposure model.
    ta_model : ImageModel or None
        Target acquisition exposure model. If None, only the science exposure is included.
    asn_fname : str, optional
        Filename for the association file.

    Returns
    -------
    str
        Path to the created association file.
    """
    sci_fname = "science.fits"
    sci_model.save(sci_fname)
    members = [{"expname": sci_fname, "exptype": "science"}]

    if ta_model is not None:
        ta_fname = "ta_image.fits"
        ta_model.save(ta_fname)
        members.append({"expname": ta_fname, "exptype": "target_acquisition"})

    asn = {
        "asn_type": "spec2",
        "asn_id": "test_ta",
        "asn_pool": "pool_id",
        "products": [
            {
                "name": "test_product",
                "members": members,
            }
        ],
    }
    with Path(asn_fname).open("w") as f:
        json.dump(asn, f)
    return asn_fname
