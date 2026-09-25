"""Mock spectral data in cal format."""

import numpy as np
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.assign_wcs.tests.test_miri import (
    create_datamodel_cube,
    create_hdul,
    create_hdul_lrs_slitless,
)
from jwst.assign_wcs.tests.test_nirspec import (
    create_nirspec_fs_file,
    create_nirspec_ifu_file,
    create_nirspec_mos_file,
)
from jwst.extract_2d.extract_2d_step import Extract2dStep

__all__ = [
    "miri_lrs_slit_cal_model",
    "miri_lrs_slitless_cal_model",
    "miri_mrs_cal_model",
    "nirspec_ifu_cal_model",
    "nirspec_mos_cal_model",
    "nirspec_slit_cal_model",
]


def miri_lrs_slit_cal_model():
    """
    Create a mock MIRI LRS FS model.

    All data arrays are zero-filled.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.SlitModel`
        The LRS slit datamodel.
    """
    hdul = create_hdul("MIRIMAGE", "ANY", "ANY")
    model = datamodels.ImageModel(hdul)
    hdul.close()

    shape = (1024, 1032)
    model.data = np.zeros(shape)
    model.err = np.zeros(shape)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = np.zeros(shape)
    model.var_rnoise = np.zeros(shape)

    # Add metadata needed for LRS FS
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    model.meta.wcsinfo.v3yangle = 0.0
    model.meta.wcsinfo.vparity = -1
    model.meta.dither.x_offset = 0.0
    model.meta.dither.y_offset = 0.0

    # Assign WCS
    model = AssignWcsStep.call(model)
    model = datamodels.SlitModel(model)

    return model


def miri_lrs_slitless_cal_model():
    """
    Create a mock MIRI LRS slitless model.

    All data arrays are zero-filled.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.SlitModel`
        The LRS slitless datamodel.
    """
    shape = (5, 416, 72)
    hdul = create_hdul_lrs_slitless()
    cube_model = create_datamodel_cube(hdul, shape)
    hdul.close()

    model = datamodels.SlitModel(cube_model)
    cube_model.close()

    model.data = np.zeros(shape)
    model.err = np.zeros(shape)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = np.zeros(shape)
    model.var_rnoise = np.zeros(shape)

    return model


def miri_mrs_cal_model(detector="MIRIFUSHORT", channel="12", band="SHORT", shape=(1024, 1032)):
    """
    Create a mock MIRI MRS model.

    Data, error, variance, and DQ planes are populated with flat data values.
    The data array is set to 1.0, error to 0.01, DQ and variances to 0.0.

    Parameters
    ----------
    detector : str, optional
        Detector name.
    channel : str, optional
        Channel name.
    band : str, optional
        Band name.
    shape : tuple of int
        Data shape.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS datamodel.
    """
    hdul = create_hdul(detector=detector, channel=channel, band=band)
    model = datamodels.IFUImageModel(hdul)
    hdul.close()

    # Add data before calling AssignWCS: the s_region depends on it existing
    model.data = np.ones(shape)
    model = AssignWcsStep.call(model)

    model.err = 0.01 * model.data
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = np.zeros(shape)
    model.var_rnoise = np.zeros(shape)
    model.var_flat = np.zeros(shape)
    return model


def nirspec_ifu_cal_model(wcs_style="coordinates"):
    """
    Create a mock NIRSpec IFU model.

    Data, error, variance, and DQ planes are populated with flat data values.
    Flat variance is not populated. Pathloss point is populated.

    The data array is set to 1.0, error to 0.01, DQ and variances to 0.0.

    Parameters
    ----------
    wcs_style : {"coordinates", "slice"}, optional
        WCS style to create.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The IFU datamodel.
    """
    shape = (2048, 2048)
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    model = datamodels.IFUImageModel(hdul)
    hdul.close()

    # assign a WCS
    if wcs_style == "coordinates":
        model = AssignWcsStep.call(model, nrs_ifu_slice_wcs=False)
    else:
        model = AssignWcsStep.call(model, nrs_ifu_slice_wcs=True)

    # add flat data
    model.data = np.ones(shape)
    model.err = 0.01 * model.data
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = np.zeros(shape)
    model.var_rnoise = np.zeros(shape)
    model.pathloss_point = np.zeros(shape)

    return model


def nirspec_mos_cal_model():
    """
    Create a mock NIRSpec MOS model.

    Calls assign_wcs and extract_2d.

    The data, dq, and variances are set to 0.0, error to 0.01.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The MOS datamodel.
    """
    hdul = create_nirspec_mos_file()
    model = datamodels.ImageModel(hdul)
    hdul.close()

    msaconfl = get_pkg_data_filename("data/msa_configuration.fits", package="jwst.assign_wcs.tests")
    model.meta.instrument.msa_metadata_file = msaconfl
    model.meta.instrument.msa_metadata_id = 12

    shape = (2048, 2048)
    model.data = np.zeros(shape)
    model.err = np.full(shape, 0.01)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = np.zeros(shape)
    model.var_rnoise = np.zeros(shape)
    model = AssignWcsStep.call(model)
    model = Extract2dStep.call(model)

    for slit in model.slits:
        slit.meta.photometry.pixelarea_steradians = 1.0
        slit.meta.photometry.pixelarea_arcsecsq = 1.0
        slit.meta.bunit_data = "MJy"

    return model


def nirspec_slit_cal_model():
    """
    Create a mock NIRSpec FS model.

    Calls assign_wcs and extract_2d.

    The data array is set to 1.0, error to 0.01, DQ and variances to 0.0.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The FS datamodel.
    """
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    model = datamodels.ImageModel(hdul)
    hdul.close()

    shape = (2048, 2048)
    model.data = np.ones(shape)
    model.err = model.data * 0.01
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.var_poisson = np.zeros(shape)
    model.var_rnoise = np.zeros(shape)
    model = AssignWcsStep.call(model)
    model = Extract2dStep.call(model)

    for slit in model.slits:
        slit.meta.photometry.pixelarea_steradians = 1.0
        slit.meta.photometry.pixelarea_arcsecsq = 1.0
        if slit.name == model.meta.instrument.fixed_slit:
            slit.meta.bunit_data = "MJy"
        else:
            slit.meta.bunit_data = "MJy/sr"

    return model
