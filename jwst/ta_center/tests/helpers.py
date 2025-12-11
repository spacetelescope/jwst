import json
from pathlib import Path

import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import models

from jwst.pathloss.pathloss import calculate_pathloss_vector
from jwst.ta_center.ta_center import JWST_DIAMETER, PIXSCALE
from jwst.ta_center.ta_center_step import _get_wavelength

X_REF_SLIT = 326.13
Y_REF_SLIT = 300.7
X_REF_SLITLESS = 38.5
Y_REF_SLITLESS = 829.0
MIRI_DETECTOR_SHAPE = (1024, 1032)  # (ny, nx) for MIRI imager


def make_pathloss_model():
    """
    Create a mock MIRI LRS pathloss reference file.

    Returns
    -------
    str
        Path to the file.
    """
    pathloss_model = dm.MirLrsPathlossModel()

    # Define wavelength grid over MIRI LRS wavelength range
    n_wavelengths = 50
    wavelengths = np.linspace(5.0, 14.0, n_wavelengths).astype(np.float32)

    # Define spatial grids - slit is long in x, narrow in y
    n_x = 51
    n_y = 21
    y_grid = np.linspace(-0.5, 0.5, n_y)  # arcsec

    # Define pathloss data. This will be a slight Gaussian drop-off from the center
    # until +/- 5 pixels in y, then it'll go to zero. Flat in x (along the slit).
    pathloss_data = np.zeros((n_wavelengths, n_y, n_x), dtype=np.float32)
    y_offsets = y_grid[:, np.newaxis]  # Shape: (n_y, 1), in arcsec
    pathloss_2d = np.exp(-2 * y_offsets**2)  # Gaussian in y, broadcast to x

    # apply the cutoff in y direction
    y_pixel_scale = (y_grid[-1] - y_grid[0]) / (n_y - 1)
    cutoff_arcsec = 5 * y_pixel_scale  # ~0.25 arcsec
    pathloss_2d[np.abs(y_offsets) > cutoff_arcsec] = 0.0
    pathloss_data[:] = pathloss_2d

    # Create the pathloss table as a structured numpy array
    # The schema requires columns: wavelength (1D), pathloss (2D), pathloss_err (2D)
    dtype = [
        ("wavelength", np.float32),
        ("pathloss", np.float32, (n_y, n_x)),
        ("pathloss_err", np.float32, (n_y, n_x)),
    ]
    pathloss_table = np.zeros(n_wavelengths, dtype=dtype)
    pathloss_table["wavelength"] = wavelengths
    pathloss_table["pathloss"] = pathloss_data
    pathloss_table["pathloss_err"] = np.ones_like(pathloss_data) * 0.01
    pathloss_model.pathloss_table = pathloss_table

    # Set up wcsinfo
    pathloss_model.meta.wcsinfo.crpix1 = (n_x + 1) / 2.0  # Reference pixel in x
    pathloss_model.meta.wcsinfo.crpix2 = (n_y + 1) / 2.0  # Reference pixel in y
    pathloss_model.meta.wcsinfo.crval1 = 0.0  # Reference value in x (arcsec)
    pathloss_model.meta.wcsinfo.crval2 = 0.0  # Reference value in y (arcsec)

    pathloss_model.meta.wcsinfo.cdelt1 = 1  # x_pixel_scale
    pathloss_model.meta.wcsinfo.cdelt2 = 1  # y_pixel_scale

    return pathloss_model


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
    x_center = X_REF_SLITLESS + offset[0]
    y_center = Y_REF_SLITLESS + offset[1]

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
    wavelength = _get_wavelength("F1500W")

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slit reference position
    x_center = X_REF_SLIT + offset[0]
    y_center = Y_REF_SLIT + offset[1]

    # Create Airy disk model on the full detector
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data_full = airy(x, y)

    # Retrieve the pathloss model
    pathloss_model = make_pathloss_model()
    pathloss_table = pathloss_model.pathloss_table
    pathloss_wcs = pathloss_model.meta.wcsinfo

    # Only calculate pathloss in a region around the slit reference point
    slit_half_height = 15  # pixels in y direction (cross-dispersion)
    slit_half_width = 60  # pixels in x direction (along slit)
    y_min = max(0, int(Y_REF_SLIT - slit_half_height))
    y_max = min(MIRI_DETECTOR_SHAPE[0], int(Y_REF_SLIT + slit_half_height))
    x_min = max(0, int(X_REF_SLIT - slit_half_width))
    x_max = min(MIRI_DETECTOR_SHAPE[1], int(X_REF_SLIT + slit_half_width))

    # Calculate pathloss only in the slit region.
    # calculate_pathloss_vector is unfortunately not vectorized.
    slit_weights = np.zeros_like(data_full)
    for i in range(y_min, y_max):
        for j in range(x_min, x_max):
            # Calculate offset from reference center in arcsec
            dx_arcsec = (j - X_REF_SLIT) * PIXSCALE
            dy_arcsec = (i - Y_REF_SLIT) * PIXSCALE

            _, pathloss_vector, _ = calculate_pathloss_vector(
                pathloss_table["pathloss"], pathloss_wcs, dx_arcsec, dy_arcsec, calc_wave=False
            )
            # Find the wavelength index closest to our wavelength
            wave_idx = np.argmin(np.abs(pathloss_table["wavelength"] - wavelength))
            slit_weights[i, j] = pathloss_vector[wave_idx]

    # Add noise
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data_full.shape)
    data_full += noise

    # Apply slit weights to the data (multiply by pathloss correction)
    data = data_full * slit_weights

    # Add bad values near the slit source position to test their handling
    data[int(Y_REF_SLIT) + 4, int(X_REF_SLIT) + 2] = np.nan
    data[int(Y_REF_SLIT) - 3, int(X_REF_SLIT) - 4] = np.inf

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
