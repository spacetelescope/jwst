import numpy as np
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


def _nirspec_rate_base(
    detector="NRS1", exptype="NRS_MSASPEC", filter_name="F170LP", grating="G235M", lamp="NONE"
):
    model = datamodels.ImageModel()

    model.meta.instrument.name = "NIRSPEC"
    model.meta.instrument.detector = detector
    model.meta.instrument.filter = filter_name
    model.meta.instrument.grating = grating
    model.meta.instrument.lamp_state = lamp

    model.meta.observation.time = "8:59:37"
    model.meta.observation.date = "2016-09-05"
    model.meta.observation.program_number = "1234"

    model.meta.exposure.type = exptype

    model.meta.dither.position_number = 1

    model.meta.wcsinfo.instance.update(DEFAULT_WCS_KW)

    return model


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
    model = _nirspec_rate_base(filter_name=filter_name, exptype=exptype, grating="MIRROR")
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
    model = _nirspec_rate_base(filter_name=filter_name, grating=grating, exptype="NRS_MSASPEC")
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
    model = _nirspec_rate_base(
        filter_name=filter_name, grating=grating, lamp=lamp, exptype="NRS_FIXEDSLIT"
    )

    # Add some more metadata needed for FS
    model.meta.wcsinfo.crval3 = 0
    model.meta.wcsinfo.wcsaxes = 3
    model.meta.wcsinfo.ctype3 = "WAVE"
    model.meta.instrument.gwa_xtilt = 0.3316612243652344
    model.meta.instrument.gwa_ytilt = 0.1260581910610199
    model.meta.instrument.fixed_slit = "S200A1"
    model.meta.subarray.name = "FULL"

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
    base_model = _nirspec_rate_base(
        detector=detector, filter_name=filter_name, grating=grating, exptype="NRS_IFU"
    )
    model = datamodels.IFUImageModel(base_model)

    # Add some more metadata needed for IFU
    model.meta.observation.date = "2026-01-01"  # chromcorr CRDS selector requires date > launch
    model.meta.wcsinfo.crval3 = 0
    model.meta.wcsinfo.wcsaxes = 3
    model.meta.wcsinfo.ctype3 = "WAVE"
    model.meta.instrument.gwa_xtilt = gwa_xtil
    model.meta.instrument.gwa_ytilt = gwa_ytil
    if gwa_tilt is not None:
        model.meta.instrument.gwa_tilt = gwa_tilt

    model.data = np.zeros((10, 10), dtype=np.float32)

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = nirspec.create_pipeline(model, ref, SLIT_Y_RANGE)
        model.meta.wcs = wcs.WCS(pipeline)

        slits = list(range(30))
        model.meta.wcs.bounding_box = nirspec.generate_compound_bbox(model, slits)

    return model
