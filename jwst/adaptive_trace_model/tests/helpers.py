import warnings

import numpy as np

from jwst.tests import spec_cal_helpers

__all__ = [
    "miri_lrs_slit_model_with_source",
    "miri_lrs_slitless_model_with_source",
    "miri_mrs_model_with_source",
    "nirspec_ifu_model_with_source",
    "nirspec_mos_model_with_source",
    "nirspec_slit_model_with_source",
    "nirspec_slit_model_with_source_and_nod",
    "profile_1d",
]


def miri_lrs_slit_model_with_source():
    """
    Create a mock MIRI LRS FS model with a simple spectral source in the data array.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.SlitModel`
        The LRS slit datamodel.
    """
    model = spec_cal_helpers.miri_lrs_slit_cal_model()
    model.data *= np.nan

    shape = model.data.shape
    ysize, xsize = shape[-2:]
    x, y = np.meshgrid(np.arange(xsize), np.arange(ysize))
    _, _, lam = model.meta.wcs(x, y)

    region_map = (~np.isnan(lam)).astype(int)
    _add_source(model, region_map, along_x=False)
    model.err[:] = 0.01

    return model


def miri_lrs_slitless_model_with_source():
    """
    Create a mock MIRI LRS slitless model with a simple spectral source in the data array.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.SlitModel`
        The LRS slitless datamodel.
    """
    model = spec_cal_helpers.miri_lrs_slitless_cal_model()
    model.data *= np.nan

    shape = model.data.shape
    ysize, xsize = shape[-2:]
    x, y = np.meshgrid(np.arange(xsize), np.arange(ysize))
    _, _, lam = model.meta.wcs(x, y)

    region_map = (~np.isnan(lam)).astype(int)
    _add_source(model, region_map, along_x=False)

    return model


def miri_mrs_model_with_source():
    """
    Create a mock MIRI MRS model with a simple spectral source in the data array.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The MRS datamodel.
    """
    model = spec_cal_helpers.miri_mrs_cal_model()
    model.data *= np.nan

    # add a simple source to each slice
    det2ab_transform = model.meta.wcs.get_transform("detector", "alpha_beta")
    region_map = det2ab_transform.label_mapper.mapper
    _add_source(model, region_map, along_x=False, bright_factor=10)

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
    model = spec_cal_helpers.nirspec_ifu_cal_model(wcs_style=wcs_style)

    # add a simple source
    shape = model.data.shape
    model.data = np.full(shape, np.nan)
    region_map = model.regions
    _add_source(model, region_map, along_x=True, bright_factor=1000)
    model.err = 0.01 * model.data

    return model


def nirspec_mos_model_with_source():
    """
    Create a mock NIRSpec MOS model with a simple spectral source in the data array.

    Calls assign_wcs and extract_2d.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The MOS datamodel.
    """
    model = spec_cal_helpers.nirspec_mos_cal_model()

    for slit in model.slits:
        region_map = (~np.isnan(slit.wavelength)).astype(int)
        _add_source(slit, region_map)

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
    model = spec_cal_helpers.nirspec_slit_cal_model()
    for slit in model.slits:
        region_map = (~np.isnan(slit.wavelength)).astype(int)
        slit.data *= 0
        _add_source(slit, region_map)

    return model


def nirspec_slit_model_with_source_and_nod():
    """
    Create a mock NIRSpec FS model with a source and a negative nod in the data array.

    Calls assign_wcs and extract_2d.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The FS datamodel.
    """
    model = spec_cal_helpers.nirspec_slit_cal_model()
    for slit in model.slits:
        region_map = (~np.isnan(slit.wavelength)).astype(int)
        slit.data *= 0
        slit.err[:] = 1e-5
        _add_source(
            slit,
            region_map,
            center_offset=0.3,
            amplitude=10,
            bright_factor=1.0,
            width=1.5,
            overwrite_data=False,
        )
        _add_source(
            slit,
            region_map,
            center_offset=-0.3,
            amplitude=-10,
            bright_factor=1.0,
            width=1.5,
            overwrite_data=False,
        )

    return model


def profile_1d(xvec, center_offset=0.0, amplitude=0.1, baseline=1.0, width=2.0, along_x=True):
    """
    Make a smooth 1D Gaussian profile.

    Parameters
    ----------
    xvec : ndarray
        X-values for the profile.
    center_offset : float, optional
        Offset for centering the Gaussian, given as a fraction of the mean ``xvec`` value.
    amplitude : float, optional
        Amplitude for the Gaussian.
    baseline : float, optional
        Background level to add.
    width : float, optional
        Gaussian width in pixels.

    Returns
    -------
    yvec : ndarray
        Gaussian y-values for the profile, centered on the middle of the ``xvec`` array.
    """
    if xvec.ndim == 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            if along_x:
                center = np.nanmean(xvec, axis=0)[None, :]
            else:
                center = np.nanmean(xvec, axis=1)[:, None]
    else:
        center = np.nanmean(xvec)
    center += center_offset * center
    peak = amplitude * np.exp(-0.5 * ((xvec - center) / width) ** 2)
    return peak + baseline


def _add_source(
    model,
    region_map,
    along_x=True,
    bright_factor=10.0,
    overwrite_data=True,
    center_offset=0.0,
    amplitude=1.0,
    baseline=0.0,
    width=2.0,
):
    ysize, xsize = model.data.shape[-2:]
    x, y = np.meshgrid(np.arange(xsize), np.arange(ysize))
    slice_numbers = np.unique(region_map[region_map > 0])
    for slice_num in slice_numbers:
        indx = region_map == slice_num
        if along_x:
            slice_y = y.astype(float)
            slice_y[~indx] = np.nan
            source = profile_1d(
                slice_y,
                center_offset=center_offset,
                amplitude=amplitude,
                baseline=baseline,
                width=width,
                along_x=along_x,
            )
        else:
            slice_x = x.astype(float)
            slice_x[~indx] = np.nan
            source = profile_1d(
                slice_x,
                center_offset=center_offset,
                amplitude=amplitude,
                baseline=baseline,
                width=width,
                along_x=along_x,
            )
        if overwrite_data:
            model.data[..., indx] = source[indx]
        else:
            model.data[..., indx] += source[indx]

        # Make one slice brighter, for threshold tests
        if slice_num == slice_numbers[len(slice_numbers) // 2]:
            model.data[..., indx] *= bright_factor
