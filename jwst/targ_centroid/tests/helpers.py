import json
from pathlib import Path

import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import models

from jwst.assign_wcs import AssignWcsStep

# Constants for MIRI imager and LRS slit dimensions
PIXSCALE = 0.11  # arcsec/pixel for MIRI imager
JWST_DIAMETER = 6.5  # meters
SLIT_REF = (326.13, 300.7)  # for fixed-slit data
FULL_SHAPE = (1024, 1032)  # (ny, nx) for MIRI imager full subarray
SLITLESS_REF = (38.5, 301.0)  # for slitless data
SLITLESSPRISM_SHAPE = (72, 416)  # for SLITLESSPRISM subarray


__all__ = [
    "make_slit_data",
    "make_slitless_data",
    "make_empty_lrs_model",
    "make_ta_model",
    "make_ta_association",
]


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
    # Create coordinate grids using slitless subarray size
    y, x = np.mgrid[0 : SLITLESSPRISM_SHAPE[1], 0 : SLITLESSPRISM_SHAPE[0]]

    # Calculate diffraction-limited Airy disk radius
    # First zero of Airy disk: theta = 1.22 * lambda / D (in radians)
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slitless reference position
    x_center = SLITLESS_REF[0] + offset[0]
    y_center = SLITLESS_REF[1] + offset[1]

    # Simulate source as an Airy disk model at diffraction-limited resolution
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data = airy(x, y)

    # Add some noise with fixed seed for reproducibility
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data.shape)
    data += noise

    # Add NaN values near the slitless source position to test that this still works ok
    data[int(SLITLESS_REF[1]) + 5, int(SLITLESS_REF[0]) + 3] = np.nan
    data[int(SLITLESS_REF[1]) - 4, int(SLITLESS_REF[0]) - 2] = np.inf

    return data


def _add_metadata(model, exptype):
    """
    Add metadata to science or TA confirm data.

    Parameters
    ----------
    model : ImageModel
        Input data model.
    exptype : str
        Exposure type of the associated science data.
    """
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "F1500W"
    model.meta.exposure.type = exptype

    model.meta.observation.date = "2025-11-10"
    model.meta.observation.time = "12:00:00"

    # Add subarray metadata
    if exptype != "MIR_LRS-SLITLESS":
        model.meta.subarray.name = "FULL"
        model.meta.subarray.xstart = 1  # 1-indexed FITS convention
        model.meta.subarray.ystart = 1
        model.meta.subarray.xsize = FULL_SHAPE[1]  # nx
        model.meta.subarray.ysize = FULL_SHAPE[0]  # ny
        model.meta.wcsinfo.v2_ref = -414.848
        model.meta.wcsinfo.v3_ref = -400.426
    else:
        # slitless data are only valid in a subarray, otherwise assign_wcs errors out
        model.meta.subarray.name = "SLITLESSPRISM"
        model.meta.subarray.xstart = 1  # 1-indexed FITS convention
        model.meta.subarray.ystart = 529
        model.meta.subarray.xsize = SLITLESSPRISM_SHAPE[1]  # nx
        model.meta.subarray.ysize = SLITLESSPRISM_SHAPE[0]  # ny
        model.meta.wcsinfo.v2_ref = -378.600
        model.meta.wcsinfo.v3_ref = -344.752

    # add wcsinfo
    model.meta.wcsinfo.roll_ref = 92.69243893875218
    model.meta.wcsinfo.ra_ref = 108.63038309653392
    model.meta.wcsinfo.dec_ref = 13.860413558600543
    model.meta.wcsinfo.v3yangle = 4.75797
    model.meta.wcsinfo.vparity = -1

    # dither info
    model.meta.dither.x_offset = 0.0
    model.meta.dither.y_offset = 0.0


def make_ta_model(data, exptype):
    """
    Create a MIRI TA verification image model.

    Parameters
    ----------
    data : ndarray
        The image data array.
    exptype : str
        Exposure type of the associated science data

    Returns
    -------
    ImageModel
        MIRI TA verification image model with metadata.
    """
    model = dm.ImageModel()
    model.data = data
    _add_metadata(model, exptype)

    # override exptype as taconfirm
    model.meta.exposure.type = "MIR_TACONFIRM"

    return model


def make_slit_data(offset):
    """
    Make a fake MIRI LRS slit TA verification image.

    The output data are zero everywhere except in the slit region.
    An Airy disk PSF is created at the source position plus offset,
    noise is added, and then the slit mask is applied.

    Parameters
    ----------
    offset : tuple
        (x,y) offset of the source from the reference center in pixels.

    Returns
    -------
    model : ImageModel
        Simulated slit TA image model.
    """
    wavelength = 15.0

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slit reference position
    x_center = SLIT_REF[0] + offset[0]
    y_center = SLIT_REF[1] + offset[1]

    # Create Airy disk model on the full detector
    y, x = np.mgrid[0 : FULL_SHAPE[0], 0 : FULL_SHAPE[1]]
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data_full = airy(x, y)

    # Make a slit mask
    slit_half_height = 8  # pixels in y direction (cross-dispersion) very approximate
    slit_half_width = 60  # pixels in x direction (along slit)
    y_min = max(0, int(SLIT_REF[1] - slit_half_height))
    y_max = min(FULL_SHAPE[0], int(SLIT_REF[1] + slit_half_height))
    x_min = max(0, int(SLIT_REF[0] - slit_half_width))
    x_max = min(FULL_SHAPE[1], int(SLIT_REF[0] + slit_half_width))
    slit_mask = np.zeros(FULL_SHAPE, dtype=bool)
    slit_mask[y_min:y_max, x_min:x_max] = True

    # Add noise
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data_full.shape)
    data_full += noise

    # Apply slit mask to the data
    data = data_full * slit_mask

    # Add bad values near the slit source position to test their handling
    data[int(SLIT_REF[1]) + 4, int(SLIT_REF[0]) + 2] = np.nan
    data[int(SLIT_REF[1]) - 3, int(SLIT_REF[0]) - 4] = np.inf

    model = make_ta_model(data, "MIR_LRS-FIXEDSLIT")
    return model


def make_empty_lrs_model(exptype):
    """
    Create empty SlitModel with metadata shared between slit and slitless modes for LRS.

    Parameters
    ----------
    exptype : str
        Exposure type, one of MIR_LRS-FIXEDSLIT or MIR_LRS-SLITLESS

    Returns
    -------
    SlitModel
        MIRI LRS model with required metadata.
    """
    model = dm.ImageModel((10, 10))
    _add_metadata(model, exptype)

    # Assign a WCS, then convert to SlitModel as would be done by calwebb_spec2
    # must be done in this order because AssignWcsStep doesn't like SlitModels
    model = AssignWcsStep.call(model)
    return dm.SlitModel(model)


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
