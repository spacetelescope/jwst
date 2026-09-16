import numpy as np
from stdatamodels.jwst.datamodels import MaskModel

from jwst.lib.tests.helpers import make_rawramp

__all__ = [
    "make_rampmodel",
    "make_dq_arrays",
]


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
    dm_ramp = make_rawramp("MIRI", nints, ngroups, ysize, xsize, 1, 1)

    # Add the detector
    dm_ramp.meta.instrument.detector = "MIRIMAGE"

    return dm_ramp


def make_dq_arrays(ysize, xsize):
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
    # create a DQ array to use in a mask model
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
