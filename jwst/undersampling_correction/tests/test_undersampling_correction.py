import numpy as np

from jwst.datamodels import RampModel
from jwst.datamodels import dqflags

from jwst.undersampling_correction.undersampling_correction import undersampling_correction

import numpy.testing as npt

test_dq_flags = dqflags.pixel
GOOD = test_dq_flags["GOOD"]
DNU = test_dq_flags["DO_NOT_USE"]


def test_pix_0():
    """
    Having all data in ramp below the signal threshold, the only DNU
    in the output GROUPDQ should be those propagated from the input.
    """
    ngroups, nints, nrows, ncols = set_scalars()
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    signal_threshold = 30000.

    # Populate pixel-specific SCI and GROUPDQ arrays
    # Set all SCI to be below the signal threshold, and some input
    # GROUPDQ to be DNU
    ramp_model.data[0, :, 0, 0] = np.array((signal_threshold - 100.), dtype=np.float32)
    ramp_model.groupdq[0, :, 0, 0] = [GOOD, DNU, GOOD, DNU, DNU, GOOD, GOOD, GOOD, GOOD, DNU]
    in_gdq = ramp_model.groupdq

    out_model = undersampling_correction(ramp_model, signal_threshold)
    out_gdq = out_model.groupdq

    npt.assert_array_equal(out_gdq, in_gdq)


def test_pix_1():
    """
    All input GROUPDQ = 'GOOD', Some ramp data exceed the signal threshold, so the only DNU
    in the output GROUPDQ should be the groups exceeding the signal threshold,
    """
    ngroups, nints, nrows, ncols = set_scalars()
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    signal_threshold = 20000.

    # Populate pixel-specific SCI and GROUPDQ arrays
    # Set seom groups' SCI to be above the signal threshold, and all input
    # GROUPDQ to be GOOD
    ramp_model.data[0, 1, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.data[0, 3, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.data[0, 8, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.groupdq[0, :, 0, 0] = [GOOD] * ngroups

    true_out_gdq = ramp_model.groupdq.copy()
    true_out_gdq[0, :, 0, 0] = [GOOD] * ngroups
    # True output DNU are at exceedance groups
    true_out_gdq[0, 1, 0, 0] = DNU
    true_out_gdq[0, 3, 0, 0] = DNU
    true_out_gdq[0, 8, 0, 0] = DNU

    out_model = undersampling_correction(ramp_model, signal_threshold)
    out_gdq = out_model.groupdq

    npt.assert_array_equal(out_gdq, true_out_gdq)


def set_scalars():
    """
    Set needed scalars for the size of the dataset,
    """
    ngroups = 10
    nints = 1
    nrows = 1
    ncols = 1

    return ngroups, nints, nrows, ncols


def create_mod_arrays(ngroups, nints, nrows, ncols):
    """
    For an input datacube (NIRISS), create arrays having
    the specified dimensions for the pixel DQ, the group DQ, and the
    ERR extensions, and create datamodels for the ramp
    """
    err = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint8)

    # Create and populate ramp model
    ramp_model = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq)
    ramp_model.meta.instrument.name = 'NIRISS'
    ramp_model.meta.instrument.detector = 'NRS1'

    ramp_model.meta.subarray.name = 'FULL'
    ramp_model.meta.subarray.xstart = 1
    ramp_model.meta.subarray.ystart = 1
    ramp_model.meta.subarray.xsize = ncols
    ramp_model.meta.subarray.ysize = nrows

    ramp_model.meta.exposure.ngroups = ngroups
    ramp_model.meta.exposure.nframes = 1
    ramp_model.meta.exposure.groupgap = 0
    ramp_model.meta.exposure.drop_frames1 = 0
    ramp_model.meta.exposure.type = 'NIS_IMAGE'

    return ramp_model, pixdq, gdq, err
