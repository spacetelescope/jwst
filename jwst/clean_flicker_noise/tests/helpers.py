import numpy as np
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.tests.test_nirspec import create_nirspec_fs_file, create_nirspec_ifu_file
from jwst.msaflagopen.tests.test_msa_open import make_nirspec_mos_model

__all__ = [
    "add_metadata",
    "make_small_ramp_model",
    "make_small_rate_model",
    "make_small_rateints_model",
    "make_flat_model",
    "make_nirspec_ifu_model",
    "make_nirspec_mos_fs_model",
    "make_nirspec_fs_model",
    "make_nirspec_mos_model",
]


def add_metadata(model, shape):
    """
    Add basic MIRI image metadata to a model.

    Parameters
    ----------
    model : DataModel
        Model to update
    shape : tuple of int
        Ramp data shape (nints, ngroups, ny, nx).
    """
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "F480M"
    model.meta.observation.date = "2015-10-13"
    model.meta.observation.time = "00:00:00"
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.exposure.group_time = 1.0
    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = shape[3]
    model.meta.subarray.ysize = shape[2]
    model.meta.exposure.frame_time = 1.0
    model.meta.exposure.ngroups = shape[1]
    model.meta.exposure.group_time = 1.0
    model.meta.exposure.nints = shape[0]
    model.meta.exposure.nframes = 1
    model.meta.exposure.groupgap = 0
    model.meta.exposure.readpatt = "FASTR1"
    model.meta.subarray.slowaxis = 2


def make_small_ramp_model(shape=(3, 5, 10, 10)):
    """
    Make a small ramp model with MIRI image metadata.

    Parameters
    ----------
    shape : tuple of int, optional
        Ramp data shape (nints, ngroups, ny, nx).

    Returns
    -------
    RampModel
        A ramp model with specified shape.
    """
    rampmodel = datamodels.RampModel(shape)
    add_metadata(rampmodel, shape)

    # Make data with a constant rate
    for group in range(shape[1]):
        rampmodel.data[:, group, :, :] = group

    return rampmodel


def make_small_rate_model(shape=(3, 5, 10, 10)):
    """
    Make a small rate model with MIRI image metadata.

    Parameters
    ----------
    shape : tuple of int, optional
        Ramp data shape (nints, ngroups, ny, nx). The nints and
        ngroups values are used to update the metadata. The
        image will have shape (ny, nx).

    Returns
    -------
    ImageModel
        An rate model with specified shape.
    """
    ratemodel = datamodels.ImageModel(shape[2:])
    add_metadata(ratemodel, shape)
    ratemodel.data[:] = 1.0
    return ratemodel


def make_small_rateints_model(shape=(3, 5, 10, 10)):
    """
    Make a small rateints model with MIRI image metadata.

    Parameters
    ----------
    shape : tuple of int, optional
        Ramp data shape (nints, ngroups, ny, nx). The ngroups
        value is used to update the metadata. The cube will have
        shape (nints, ny, nx).

    Returns
    -------
    CubeModel
        A rateints model with specified shape.
    """
    ratemodel = datamodels.CubeModel((shape[0], shape[2], shape[3]))
    add_metadata(ratemodel, shape)
    ratemodel.data[:] = 1.0
    return ratemodel


def make_flat_model(model, shape=(10, 10), value=None):
    """
    Make a flat model to match a science model.

    Parameters
    ----------
    model : RampModel, ImageModel, or CubeModel
        The science model to match.
    shape : tuple of int, optional
        Must be either full size for the instrument or else match
        the input model.
    value : float or None, optional
        If provided, the flat array is filled with this value. If
        None, the array is filled with 2D gradient values (via numpy arange).

    Returns
    -------
    FlatModel
        The flat model.
    """
    # make a flat model with appropriate size and metadata
    flat = datamodels.FlatModel()
    if value is None:
        flat.data = np.arange(shape[0] * shape[1], dtype=float).reshape(shape)
    else:
        flat.data = np.full(shape, value)

    # add required metadata
    flat.meta.description = "test"
    flat.meta.reftype = "test"
    flat.meta.author = "test"
    flat.meta.pedigree = "test"
    flat.meta.useafter = "test"

    # copy any other matching metadata
    flat.update(model)

    # make sure shape keys match input
    flat.meta.subarray.xsize = shape[1]
    flat.meta.subarray.ysize = shape[0]

    return flat


def make_nirspec_ifu_model(shape=(2048, 2048)):
    """
    Make a NIRSpec IFU model with realistic metadata.

    Metadata is sufficient to allow a WCS to be assigned.

    Parameters
    ----------
    shape : tuple of int, optional
        Output data shape. Must be (2048, 2048) for ``assign_wcs``
        to succeed.

    Returns
    -------
    IFUImageModel
        The IFU image model.
    """
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    hdul["SCI"].data = np.ones(shape, dtype=float)
    rate_model = datamodels.IFUImageModel(hdul)
    hdul.close()
    return rate_model


def make_nirspec_mos_fs_model():
    """
    Make a NIRSpec MOS/FS model with realistic metadata.

    Metadata is sufficient to allow a WCS to be assigned.

    Returns
    -------
    ImageModel
        The MOS/FS image model.
    """
    mos_model = make_nirspec_mos_model()
    mos_model.meta.instrument.msa_metadata_file = get_pkg_data_filename(
        "data/msa_fs_configuration.fits", package="jwst.assign_wcs.tests"
    )
    return mos_model


def make_nirspec_fs_model():
    """
    Make a NIRSpec FS model with realistic metadata.

    Metadata is sufficient to allow a WCS to be assigned.

    Returns
    -------
    ImageModel
        The MOS/FS image model.
    """
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    hdul["SCI"].data = np.ones((2048, 2048), dtype=float)
    rate_model = datamodels.ImageModel(hdul)
    hdul.close()

    # add the slow axis and subarray information
    rate_model.meta.subarray.slowaxis = 1
    rate_model.meta.subarray.fastaxis = 2
    rate_model.meta.subarray.xstart = 1
    rate_model.meta.subarray.ystart = 1
    rate_model.meta.subarray.xsize = 2048
    rate_model.meta.subarray.ysize = 2048
    return rate_model


def make_niriss_rate_model(shape=None):
    """
    Make a NIRISS image rate model with basic metadata.

    Parameters
    ----------
    shape : tuple of int or None, optional
        If not provided, the default shape in `make_small_rate_model`
        is used.

    Returns
    -------
    ImageModel
        A NIRISS image model.
    """
    if shape is not None:
        model = make_small_rate_model(shape=shape)
    else:
        model = make_small_rate_model()
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "CLEAR"
    model.meta.exposure.type = "NIS_IMAGE"
    model.meta.subarray.slowaxis = -1
    return model


def make_nircam_rate_model(shape=None):
    """
    Make a NIRCam image rate model with basic metadata.

    Parameters
    ----------
    shape : tuple of int or None, optional
        If not provided, the default shape in `make_small_rate_model`
        is used.

    Returns
    -------
    ImageModel
        A NIRCam image model.
    """
    if shape is not None:
        model = make_small_rate_model(shape=shape)
    else:
        model = make_small_rate_model()
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.detector = "NRCBLONG"
    model.meta.instrument.filter = "F300M"
    model.meta.instrument.pupil = "CLEAR"
    model.meta.exposure.type = "NRC_IMAGE"
    model.meta.subarray.slowaxis = -2
    return model


def make_nrs_fs_full_ramp():
    shape = (1, 2, 2048, 2048)
    model = datamodels.RampModel(shape)

    # Make data with a constant rate
    for group in range(shape[1]):
        model.data[:, group, :, :] = group

    # Add NIRSpec metadata from a rate model
    image_model = make_nirspec_fs_model()
    model.update(image_model)
    image_model.close()

    # Add ramp information
    model.meta.exposure.nints = 1
    model.meta.exposure.ngroups = 2
    model.meta.exposure.nframes = 1
    model.meta.exposure.groupgap = 0
    model.meta.exposure.group_time = 1.0
    model.meta.exposure.frame_time = 14.5
    model.meta.exposure.readpatt = "NRS"

    return model
