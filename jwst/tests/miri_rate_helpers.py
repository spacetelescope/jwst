import numpy as np
from gwcs import wcs
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import miri
from jwst.tests.wcs_helpers import get_reference_files

__all__ = [
    "miri_image_rate_model",
    "miri_mrs_rate_model",
    "miri_wfss_rate_model",
    "miri_lrs_slitless_rate_model",
    "miri_lrs_slit_rate_model",
]


# Default WCS information
DEFAULT_WCS_KW = {
    "wcsaxes": 3,
    "ra_ref": 165,
    "dec_ref": 54,
    "v2_ref": -8.3942412,
    "v3_ref": -5.3123744,
    "roll_ref": 37,
    "crpix1": 1024,
    "crpix2": 1024,
    "crpix3": 0,
    "cdelt1": 0.08,
    "cdelt2": 0.08,
    "cdelt3": 1,
    "ctype1": "RA---TAN",
    "ctype2": "DEC--TAN",
    "ctype3": "WAVE",
    "pc1_1": 1,
    "pc1_2": 0,
    "pc1_3": 0,
    "pc2_1": 0,
    "pc2_2": 1,
    "pc2_3": 0,
    "pc3_1": 0,
    "pc3_2": 0,
    "pc3_3": 1,
    "cunit1": "deg",
    "cunit2": "deg",
    "cunit3": "um",
}


def _miri_rate_base(detector="MIRIMAGE", channel="ANY", band="ANY", exptype="MIR_IMAGE"):
    model = datamodels.ImageModel()
    model.meta.telescope = "JWST"

    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = detector
    model.meta.instrument.channel = channel
    model.meta.instrument.band = band

    model.meta.observation.time = "8:59:37"
    model.meta.observation.date = "2017-09-05"

    model.meta.exposure.type = exptype

    model.meta.wcsinfo.instance.update(DEFAULT_WCS_KW)

    return model


def miri_image_rate_model(with_wcs=True):
    """
    Create a mock MIRI image rate model.

    The data array is zero-filled with shape (10, 10).

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The MIRI image datamodel.
    """
    model = _miri_rate_base(exptype="MIR_IMAGE")
    model.data = np.zeros((10, 10), dtype=np.float32)

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = miri.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def miri_mrs_rate_model(
    detector="MIRIFUSHORT", channel="12", band="SHORT", shape=(10, 10), with_wcs=True
):
    """
    Create a mock MIRI MRS rate model.

    The data array is zero-filled with the given shape.

    Parameters
    ----------
    detector : str, optional
        Detector name.
    channel : str, optional
        Channel name.
    band : str, optional
        Band name.
    shape : tuple of int, optional
        Shape of the data array to assign.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel` or \
            `~stdatamodels.jwst.datamodels.CubeModel`
        The MIRI datamodel. Image model if data shape is 2D; cube model if 3D.
    """
    base_model = _miri_rate_base(detector=detector, channel=channel, band=band, exptype="MIR_MRS")
    if len(shape) == 3:
        model = datamodels.CubeModel(base_model)
    else:
        model = datamodels.IFUImageModel(base_model)
    model.data = np.zeros(shape, dtype=np.float32)
    model.meta.filename = f"test{channel}{band}"

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = miri.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def miri_wfss_rate_model(with_wcs=True):
    """
    Create a mock MIRI WFSS rate model.

    The data array is zero-filled with shape (10,10).

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The MIRI image datamodel.
    """
    model = _miri_rate_base(exptype="MIR_WFSS")

    # Give the model a subarray and newer date to pick up appropriate reference files
    model.meta.subarray.name = "FULL"
    model.meta.observation.date = "2026-09-05"

    model.data = np.zeros((10, 10), dtype=np.float32)
    model.meta.filename = "test_miri_wfss"
    model.meta.instrument.filter = "P750L"

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = miri.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def miri_lrs_slitless_rate_model(with_wcs=True, subarray="SLITLESSPRISM"):
    """
    Create a mock MIRI LRS slitless rate model.

    The data array is zero-filled with shape (5, 416, 72).

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.
    subarray : str, optional
        Subarray name, either "SLITLESSPRISM" or "SLITLESSPRISM_IPS".

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        The MIRI LRS slitless datamodel.
    """
    base_model = _miri_rate_base(exptype="MIR_LRS-SLITLESS")
    model = datamodels.CubeModel(base_model)
    model.meta.subarray.name = subarray
    if subarray == "SLITLESSPRISM_IPS":
        model.meta.subarray.xstart = 13
        model.meta.subarray.ystart = 753
        model.meta.subarray.xsize = 52
        model.meta.subarray.ysize = 256
    else:
        model.meta.subarray.xstart = 1
        model.meta.subarray.ystart = 529
        model.meta.subarray.xsize = 72
        model.meta.subarray.ysize = 416

    model.data = np.zeros(
        (5, model.meta.subarray.ysize, model.meta.subarray.xsize), dtype=np.float32
    )

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = miri.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def miri_lrs_slit_rate_model(with_wcs=True):
    """
    Create a mock MIRI LRS fixed slit rate model.

    The data array is zero-filled with shape (1024, 1032).

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The MIRI LRS slit datamodel.
    """
    model = _miri_rate_base(exptype="MIR_LRS-FIXEDSLIT")
    model.data = np.zeros((1024, 1032), dtype=np.float32)

    # Add metadata needed for LRS FS
    model.meta.wcsinfo.v3yangle = 0.0
    model.meta.wcsinfo.vparity = -1
    model.meta.dither.x_offset = 0.0
    model.meta.dither.y_offset = 0.0

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = miri.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model
