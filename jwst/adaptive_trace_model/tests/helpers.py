import numpy as np
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.assign_wcs.tests.test_miri import create_hdul
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_fs_file, create_nirspec_ifu_file
from jwst.extract_2d.extract_2d_step import Extract2dStep

__all__ = [
    "miri_mrs_model",
    "miri_mrs_model_with_source",
    "nirspec_ifu_model",
    "nirspec_ifu_model_with_source",
    "nirspec_slit_model_with_source",
    "profile_1d",
]


def miri_mrs_model(detector="MIRIFUSHORT", channel="12", band="SHORT", shape=(1024, 1032)):
    """
    Create a mock MIRI MRS model.

    Data, error, variance, and DQ planes are populated with flat data values.

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


def miri_mrs_model_with_source():
    """
    Create a mock MIRI MRS model with a simple spectral source in the data array.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS datamodel.
    """
    model = miri_mrs_model()
    model.data *= np.nan

    # add a simple source to each slice
    det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
    region_map = det2ab_transform.label_mapper.mapper
    _add_source(model, region_map, along_x=False, bright_factor=10)

    return model


def nirspec_ifu_model(wcs_style="coordinates"):
    """
    Create a mock NIRSpec IFU model.

    Data, error, variance, and DQ planes are populated with flat data values.
    Flat variance is not populated. Pathloss point is populated.

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


def nirspec_ifu_model_with_source(wcs_style="coordinates"):
    """
    Create a mock NIRSpec IFU model with a simple spectral source in the data array.

    Parameters
    ----------
    wcs_style : {"coordinates", "slice"}, optional
        WCS style to create.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The IFU datamodel.
    """
    model = nirspec_ifu_model(wcs_style=wcs_style)

    # add a simple source
    shape = model.data.shape
    model.data = np.full(shape, np.nan)
    region_map = model.regions
    _add_source(model, region_map, along_x=True, bright_factor=1000)
    model.err = 0.01 * model.data

    return model


def nirspec_slit_model_with_source():
    """
    Create a mock NIRSpec FS model with a simple spectral source in the data array.

    Calls assign_wcs and extract_2d.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The FS datamodel.
    """
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    model = datamodels.ImageModel(hdul)
    hdul.close()

    shape = (2048, 2048)
    model.data = np.full(shape, np.nan)
    model.err = np.zeros(shape)
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
        region_map = (~np.isnan(slit.wavelength)).astype(int)
        _add_source(slit, region_map)

    return model


def profile_1d(xvec, amplitude=0.1, baseline=1.0):
    """
    Make a smooth 1D Gaussian profile.

    Parameters
    ----------
    xvec : ndarray
        X-values for the profile.
    amplitude : float, optional
        Amplitude for the Gaussian.
    baseline : float, optional
        Background level to add.

    Returns
    -------
    yvec : ndarray
        Gaussian y-values for the profile, centered on the middle of the ``xvec`` array.
    """
    xpos = np.mean(xvec)
    peak = amplitude * np.exp(-0.5 * ((xvec - xpos) / 2) ** 2)
    return peak + baseline


def _add_source(model, region_map, along_x=True, bright_factor=10.0):
    ysize, xsize = model.data.shape[-2:]
    x, y = np.meshgrid(np.arange(xsize), np.arange(ysize))
    slice_numbers = np.unique(region_map[region_map > 0])
    for slice_num in slice_numbers:
        indx = region_map == slice_num
        if along_x:
            model.data[..., indx] = profile_1d(y[indx], amplitude=1.0, baseline=0.0)
        else:
            model.data[..., indx] = profile_1d(x[indx], amplitude=1.0, baseline=0.0)

        # Make one slice brighter, for threshold tests
        if slice_num == slice_numbers[len(slice_numbers) // 2]:
            model.data[..., indx] *= bright_factor
