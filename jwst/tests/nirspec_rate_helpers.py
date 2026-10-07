import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from gwcs import wcs
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import nirspec
from jwst.tests.wcs_helpers import get_reference_files

__all__ = [
    "nirspec_image_rate_model",
    "nirspec_mos_rate_model",
    "nirspec_fs_rate_model",
    "nirspec_ifu_rate_model",
]

DEFAULT_WCS_KW = {
    "wcsaxes": 2,
    "ra_ref": 165,
    "dec_ref": 54,
    "v2_ref": -8.3942412,
    "v3_ref": -5.3123744,
    "roll_ref": 37,
    "crpix1": 1024,
    "crpix2": 1024,
    "cdelt1": 0.08,
    "cdelt2": 0.08,
    "ctype1": "RA---TAN",
    "ctype2": "DEC--TAN",
    "pc1_1": 1,
    "pc1_2": 0,
    "pc2_1": 0,
    "pc2_2": 1,
}

SLIT_Y_RANGE = [-0.5, 0.5]


def _nirspec_rate_hdul(
    detector="NRS1", exptype="NRS_MSASPEC", filter_name="F170LP", grating="G235M", lamp="NONE"
):
    """Create a fits HDUList instance."""
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header["instrume"] = "NIRSPEC"
    phdu.header["detector"] = detector
    phdu.header["time-obs"] = "8:59:37"
    phdu.header["date-obs"] = "2016-09-05"
    phdu.header["program"] = "1234"
    phdu.header["exp_type"] = exptype
    phdu.header["filter"] = filter_name
    phdu.header["grating"] = grating
    phdu.header["lamp"] = lamp
    phdu.header["PATT_NUM"] = 1

    scihdu = fits.ImageHDU()
    scihdu.header["EXTNAME"] = "SCI"
    scihdu.header.update(DEFAULT_WCS_KW)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def nirspec_image_rate_model(filter_name="F290LP", exptype="NRS_IMAGE", with_wcs=True):
    """
    Create a mock NIRSpec image rate model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    filter_name : str, optional
        Filter name.
    exptype : str, optional
        Exposure type.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The NIRSpec image datamodel.
    """
    hdul = _nirspec_rate_hdul(filter_name=filter_name, exptype=exptype, grating="MIRROR")
    model = datamodels.ImageModel(hdul)
    model.data = np.zeros((10, 10), dtype=np.float32)

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nirspec.create_pipeline(model, ref, SLIT_Y_RANGE)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nirspec_mos_rate_model(filter_name="F170LP", grating="G235M", with_wcs=True):
    """
    Create a mock NIRSpec MOS rate model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    filter_name : str, optional
        Filter name.
    grating : str, optional
        Grating name.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The NIRSpec MOS datamodel.
    """
    hdul = _nirspec_rate_hdul(filter_name=filter_name, grating=grating, exptype="NRS_MSASPEC")
    model = datamodels.ImageModel(hdul)
    model.data = np.zeros((10, 10), dtype=np.float32)

    msaconfl = get_pkg_data_filename("data/msa_configuration.fits", package="jwst.assign_wcs.tests")
    model.meta.instrument.msa_metadata_file = msaconfl
    model.meta.instrument.msa_metadata_id = 12

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nirspec.create_pipeline(model, ref, SLIT_Y_RANGE)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nirspec_fs_rate_model(filter_name="F100LP", grating="G140M", lamp="N/A", with_wcs=True):
    """
    Create a mock NIRSpec FS rate model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    filter_name : str, optional
        Filter name.
    grating : str, optional
        Grating name.
    lamp : str, optional
        Lamp name.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The NIRSpec FS datamodel.
    """
    image = _nirspec_rate_hdul(
        filter_name=filter_name, grating=grating, lamp=lamp, exptype="NRS_FIXEDSLIT"
    )

    # Add some more metadata needed for FS
    image[1].header["crval3"] = 0
    image[1].header["wcsaxes"] = 3
    image[1].header["ctype3"] = "WAVE"
    image[0].header["GWA_XTIL"] = 0.3316612243652344
    image[0].header["GWA_YTIL"] = 0.1260581910610199
    image[0].header["SUBARRAY"] = "FULL"
    image[0].header["FXD_SLIT"] = "S200A1"

    model = datamodels.ImageModel(image)
    model.data = np.zeros((10, 10), dtype=np.float32)

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nirspec.create_pipeline(model, ref, SLIT_Y_RANGE)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def nirspec_ifu_rate_model(
    detector="NRS1",
    filter_name="CLEAR",
    grating="PRISM",
    gwa_xtil=0.35986012,
    gwa_ytil=0.13448857,
    gwa_tilt=37.1,
    with_wcs=True,
):
    """
    Create a mock NIRSpec IFU rate model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    detector : str, optional
        Detector name.
    filter_name : str, optional
        Filter name.
    grating : str, optional
        Grating name.
    gwa_xtil : float, optional
        Grating wheel x-tilt.
    gwa_ytil : float, optional
        Grating wheel y-tilt.
    gwa_tilt : float, optional
        Grating wheel temperature.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The NIRSpec IFU datamodel.
    """
    image = _nirspec_rate_hdul(
        detector=detector, filter_name=filter_name, grating=grating, exptype="NRS_IFU"
    )

    # Add some more metadata needed for IFU
    image[0].header["date-obs"] = "2026-01-01"  # chromcorr CRDS selector requires date > launch
    image[1].header["crval3"] = 0
    image[1].header["wcsaxes"] = 3
    image[1].header["ctype3"] = "WAVE"
    image[0].header["GWA_XTIL"] = gwa_xtil
    image[0].header["GWA_YTIL"] = gwa_ytil
    if gwa_tilt is not None:
        image[0].header["GWA_TILT"] = gwa_tilt

    model = datamodels.IFUImageModel(image)
    model.data = np.zeros((10, 10), dtype=np.float32)

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nirspec.create_pipeline(model, ref, SLIT_Y_RANGE)
        model.meta.wcs = wcs.WCS(pipeline)

        slits = list(range(30))
        model.meta.wcs.bounding_box = nirspec.generate_compound_bbox(model, slits)

    return model
