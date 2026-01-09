import numpy as np
from stdatamodels.jwst import datamodels

__all__ = ["make_miri_ramp_model", "make_nircam_rate_model", "make_nirspec_ifu_rate_model"]


def make_miri_ramp_model(nints=1, ngroups=5, ysize=1024, xsize=1032):
    """
    Make a MIRI image ramp model with minimal metadata.

    Parameters
    ----------
    nints : int, optional
        Number of integrations.
    ngroups : int, optional
        Number of groups.
    ysize : int, optional
        Y size.
    xsize : int, optional
        X size.

    Returns
    -------
    ramp : `stdatamodels.jwst.datamodels.RampModel`
        A ramp model with only the data array and minimal metadata set.
    """
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)

    ramp = datamodels.RampModel(data=data)
    ramp.meta.filename = "test_miri.fits"
    ramp.meta.instrument.name = "MIRI"
    ramp.meta.instrument.detector = "MIRIMAGE"
    ramp.meta.exposure.frame_time = 1.0
    ramp.meta.exposure.groupgap = 0
    ramp.meta.exposure.group_time = 1.0
    ramp.meta.exposure.nframes = 1
    ramp.meta.exposure.nints = nints
    ramp.meta.exposure.ngroups = ngroups
    ramp.meta.exposure.readpatt = "FASTR1"
    ramp.meta.observation.date = "2024-01-01"
    ramp.meta.observation.time = "00:00:00"
    ramp.meta.subarray.name = "FULL"
    ramp.meta.subarray.xstart = 1
    ramp.meta.subarray.xsize = xsize
    ramp.meta.subarray.ystart = 1
    ramp.meta.subarray.ysize = ysize

    return ramp


def make_nircam_rate_model():
    """
    Make a NIRCam image rate model with minimal metadata.

    Returns
    -------
    rate : `stdatamodels.jwst.datamodels.ImageModel`
        The rate model.
    """
    image = datamodels.ImageModel((2048, 2048))
    image.data[:, :] = 1
    image.meta.instrument.name = "NIRCAM"
    image.meta.instrument.filter = "F200W"
    image.meta.instrument.pupil = "CLEAR"
    image.meta.exposure.type = "NRC_IMAGE"
    image.meta.observation.date = "2019-02-27"
    image.meta.observation.time = "13:37:18.548"
    image.meta.date = "2019-02-27T13:37:18.548"
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1

    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]

    image.meta.instrument.channel = "SHORT"
    image.meta.instrument.module = "A"
    image.meta.instrument.detector = "NRCA1"

    # bare minimum wcs info to get assign_wcs step to pass
    image.meta.wcsinfo.crpix1 = 693.5
    image.meta.wcsinfo.crpix2 = 512.5
    image.meta.wcsinfo.v2_ref = -453.37849
    image.meta.wcsinfo.v3_ref = -373.810549
    image.meta.wcsinfo.roll_ref = 272.3237653262276
    image.meta.wcsinfo.ra_ref = 80.54724018120017
    image.meta.wcsinfo.dec_ref = -69.5081101864959

    return image


def make_nirspec_ifu_rate_model():
    """
    Make a NIRSpec IFU rate model with minimal metadata.

    Returns
    -------
    rate : `stdatamodels.jwst.datamodels.IFUImageModel`
        The rate model.
    """
    image = datamodels.IFUImageModel((2048, 2048))
    image.data[:, :] = 1

    image.meta.instrument.name = "NIRSPEC"
    image.meta.instrument.detector = "NRS1"
    image.meta.instrument.filter = "CLEAR"
    image.meta.instrument.grating = "PRISM"
    image.meta.exposure.type = "NRS_IFU"
    image.meta.observation.date = "2019-02-27"
    image.meta.observation.time = "13:37:18.548"
    image.meta.date = "2019-02-27T13:37:18.548"
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1
    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]
    image.meta.instrument.gwa_tilt = 37.0610

    # bare minimum wcs info to get assign_wcs step to pass
    image.meta.wcsinfo.crpix1 = 693.5
    image.meta.wcsinfo.crpix2 = 512.5
    image.meta.wcsinfo.v2_ref = -453.37849
    image.meta.wcsinfo.v3_ref = -373.810549
    image.meta.wcsinfo.roll_ref = 272.3237653262276
    image.meta.wcsinfo.ra_ref = 80.54724018120017
    image.meta.wcsinfo.dec_ref = -69.5081101864959

    return image
