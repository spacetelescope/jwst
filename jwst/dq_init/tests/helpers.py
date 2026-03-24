import numpy as np
from stdatamodels.jwst.datamodels import GuiderRawModel, MaskModel, RampModel, dqflags

__all__ = [
    "make_rawramp",
    "make_rampmodel",
    "make_maskmodel",
    "make_superstripe_model",
    "one_stripe_mask",
    "make_superstripe_mask_model",
]


def make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart, exp_type=None):
    """
    Make a raw ramp model for any instrument.

    Parameters
    ----------
    instrument : str
        Instrument name.
    nints : int
        Number of integrations.
    ngroups : int
        Number of groups.
    ysize : int
        Y-size for the subarray.
    xsize : int
        X-size for the subarray.
    ystart : int
        Y-start value for the subarray.
    xstart : int
        X-start value for the subarray.
    exp_type : str, optional
        Exposure type.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The ramp model.
    """
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)

    # create a JWST datamodel
    if instrument == "FGS":
        dm_ramp = GuiderRawModel(data=data)
        dm_ramp.meta.exposure.type = exp_type
    else:
        dm_ramp = RampModel(data=data)

    dm_ramp.meta.subarray.xstart = xstart
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = ystart
    dm_ramp.meta.subarray.ysize = ysize

    # add some basic metadata
    dm_ramp.meta.instrument.name = instrument
    dm_ramp.meta.observation.date = "2026-01-01"
    dm_ramp.meta.observation.time = "00:00:00"
    dm_ramp.meta.exposure.nints = nints
    dm_ramp.meta.exposure.ngroups = ngroups
    dm_ramp.meta.exposure.nframes = 1

    return dm_ramp


def make_rampmodel(nints=1, ngroups=5, ysize=1024, xsize=1032):
    """
    Make a MIRI ramp model with sufficient metadata to run the full step.

    Parameters
    ----------
    nints : int, optional
        Number of integrations.
    ngroups : int, optional
        Number of groups.
    ysize : int, optional
        Y-size for the subarray.
    xsize : int, optional
        X-size for the subarray

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The ramp model.
    """
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

    dm_ramp.meta.instrument.name = "MIRI"
    dm_ramp.meta.instrument.detector = "MIRIMAGE"
    dm_ramp.meta.observation.date = "2018-01-01"
    dm_ramp.meta.observation.time = "00:00:00"
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp


def make_maskmodel(ysize, xsize):
    """
    Make dq and dq_def arrays for use in creating a DQ mask model.

    Parameters
    ----------
    ysize : int
        Array y-size.
    xsize : int
        Array x-size.

    Returns
    -------
    dq : ndarray
        The DQ image.
    dq_def : ndarray
        The corresponding DQ definition table.
    """
    # create a mask model for the dq_init step
    csize = (ysize, xsize)
    dq = np.zeros(csize, dtype=int)
    # define a dq_def extension
    mask = MaskModel()

    dqdef = [
        (0, 1, "DO_NOT_USE", "Bad Pixel do not use"),
        (1, 2, "DEAD", "Dead Pixel"),
        (2, 4, "HOT", "Hot pixel"),
        (3, 8, "UNRELIABLE_SLOPE", "Large slope variance"),
        (4, 16, "RC", "RC pixel"),
        (5, 32, "REFERENCE_PIXEL", "Reference Pixel"),
    ]

    dq_def = np.array((dqdef), dtype=mask.get_dtype("dq_def"))

    return dq, dq_def


def make_superstripe_model():
    """
    Make a NIRISS SOSS superstripe raw ramp model.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The ramp model.
    """
    nints = 3
    nstripe = 10
    ngroups = 5
    ysize = 256
    xsize = 208
    xstart = 1
    ystart = 1793
    model = make_rawramp("NIRISS", nints * nstripe, ngroups, ysize, xsize, ystart, xstart)

    # Add some more basic metadata
    model.meta.instrument.detector = "NIS"
    model.meta.subarray.name = "SUB204STRIPE_SOSS"
    model.meta.subarray.fastaxis = -2
    model.meta.subarray.slowaxis = -1
    model.meta.exposure.integration_start = 1
    model.meta.exposure.integration_end = nints * nstripe
    model.meta.exposure.nints = nints * nstripe

    # Add the superstripe metadata
    model.meta.subarray.multistripe_reads1 = 4
    model.meta.subarray.multistripe_reads2 = 204
    model.meta.subarray.multistripe_skips1 = 0
    model.meta.subarray.multistripe_skips2 = 0
    model.meta.subarray.num_superstripe = nstripe
    model.meta.subarray.repeat_stripe = 1
    model.meta.subarray.interleave_reads1 = 1
    model.meta.subarray.superstripe_step = 204

    return model


def one_stripe_mask(shape):
    """
    Make a DQ mask with a range of bad pixels for a single stripe.

    Parameters
    ----------
    shape : tuple of int
        Stripe data shape.

    Returns
    -------
    dq : ndarray
        DQ array matching the stripe with bad pixels.
    """
    dq = np.zeros(shape, dtype=np.uint32)
    dq[2, 100] = dqflags.pixel["DEAD"]
    dq[4, 100] = dqflags.pixel["HOT"]
    dq[6, 100] = dqflags.pixel["UNRELIABLE_SLOPE"]
    dq[8, 100] = dqflags.pixel["RC"]
    dq[10, 100] = dqflags.pixel["DO_NOT_USE"]
    dq[2, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["DEAD"]
    dq[4, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["HOT"]
    dq[6, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["UNRELIABLE_SLOPE"]
    dq[8, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["RC"]
    return dq


def make_superstripe_mask_model():
    """
    Create a MaskModel with bad pixels in a superstripe region.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.MaskModel`
        The mask model.
    """
    ysize = 2048
    xsize = 2048
    dq = np.zeros((ysize, xsize), dtype=np.uint32)

    # edit reference file with known bad pixel values
    substart = 1792
    subsize = 256
    nstripe = 10
    stripe = 204
    stripe_shape = (subsize, stripe)
    for i in range(nstripe):
        xstart = stripe * i + 4
        xstop = xstart + stripe
        dq[substart : substart + subsize, xstart:xstop] = one_stripe_mask(stripe_shape)

    # write mask model
    ref_data = MaskModel(dq=dq)
    ref_data.meta.instrument.name = "NIRISS"
    ref_data.meta.subarray.xstart = 1
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = 1
    ref_data.meta.subarray.ysize = ysize
    ref_data.meta.description = "Test description"
    ref_data.meta.reftype = "mask"
    ref_data.meta.author = "Test Author"
    ref_data.meta.pedigree = "test"
    ref_data.meta.useafter = "2024-01-01"

    return ref_data
