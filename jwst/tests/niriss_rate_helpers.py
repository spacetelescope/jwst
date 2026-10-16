import numpy as np
from gwcs import wcs
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import niriss
from jwst.tests.wcs_helpers import get_reference_files

__all__ = ["niriss_image_rate_model", "niriss_wfss_rate_model"]


# default wcs information
DEFAULT_WCS_KW = {
    "wcsaxes": 2,
    "ra_ref": 53.16255,
    "dec_ref": -27.791461111111,
    "v2_ref": -289.966976,
    "v3_ref": -697.723326,
    "roll_ref": 0,
    "crval1": 53.16255,
    "crval2": -27.791461111111,
    "crpix1": 1024.5,
    "crpix2": 1024.5,
    "cdelt1": 0.065398,
    "cdelt2": 0.065893,
    "ctype1": "RA---TAN",
    "ctype2": "DEC--TAN",
    "pc1_1": 1,
    "pc1_2": 0,
    "pc2_1": 0,
    "pc2_2": 1,
    "cunit1": "deg",
    "cunit2": "deg",
}


def _niriss_rate_base(detector="NIS", filter_name="CLEAR", exptype="NIS_IMAGE", pupil="F200W"):
    model = datamodels.ImageModel()
    model.meta.telescope = "JWST"
    model.meta.filename = "test+" + filter_name

    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = detector
    model.meta.instrument.filter = filter_name
    model.meta.instrument.pupil = pupil
    model.meta.instrument.filter_position = 354.2111

    model.meta.observation.time = "8:59:37"
    model.meta.observation.date = "2022-09-05"

    model.meta.exposure.type = exptype

    model.meta.wcsinfo.instance.update(DEFAULT_WCS_KW)

    model.data = np.zeros((2048, 2048))

    return model


def niriss_image_rate_model(with_wcs=True):
    """
    Create a mock NIRISS image rate model.

    The data array is zero-filled with shape (2048, 2048).

    Parameters
    ----------
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The NIRISS image datamodel.
    """
    model = _niriss_rate_base(filter_name="F200W", pupil="CLEAR")

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = niriss.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model


def niriss_wfss_rate_model(filter_name="GR150R", with_wcs=True):
    """
    Create a mock NIRISS WFSS rate model.

    The data array is zero-filled with shape (2048, 2048).

    Parameters
    ----------
    filter_name : str, optional
        Filter name.
    with_wcs : bool, optional
        If True, assign a WCS to the output model.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The NIRISS WFSS datamodel.
    """
    model = _niriss_rate_base(filter_name=filter_name, pupil="F200W", exptype="NIS_WFSS")

    if with_wcs:
        ref = get_reference_files(model)
        pipeline = niriss.create_pipeline(model, ref)
        model.meta.wcs = wcs.WCS(pipeline)

    return model
