import numpy as np

from stdatamodels.jwst.datamodels import RampModel
from stdatamodels.jwst.datamodels import dqflags

from jwst.undersampling_correction.undersampling_correction import undersampling_correction
from jwst.undersampling_correction.undersampling_correction_step import UndersamplingCorrectionStep

import numpy.testing as npt

test_dq_flags = dqflags.pixel
GOOD = test_dq_flags["GOOD"]
DNU = test_dq_flags["DO_NOT_USE"]
UNSA = test_dq_flags["UNDERSAMP"]
ADFL = test_dq_flags["AD_FLOOR"]
DROU = test_dq_flags["DROPOUT"]


def test_pix_0():
    """
    Having all data in ramp below the signal threshold, the only non-GOOD
    groups in the output GROUPDQ should be those DNU propagated from the input.
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
    All input GROUPDQ = 'GOOD'. Some ramp data exceed the signal threshold, so the
    only non-GOOD groups in the output GROUPDQ should be UNSA + DNU for the first
    group exceeding the signal threshold and all subsequent groups.
    """
    ngroups, nints, nrows, ncols = set_scalars()
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    signal_threshold = 20000.

    # Populate pixel-specific SCI and GROUPDQ arrays
    # Set seom groups' SCI to be above the signal threshold, and all input
    # GROUPDQ to be GOOD
    ramp_model.data[0, 3, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.data[0, 5, 0, 0] = np.array((signal_threshold - 200.), dtype=np.float32)
    ramp_model.data[0, 8, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.groupdq[0, :, 0, 0] = [GOOD] * ngroups

    true_out_gdq = ramp_model.groupdq.copy()  # all GOOD
    true_out_gdq[0, :, 0, 0] = [GOOD] * ngroups
    true_out_gdq[0, 3:, 0, 0] = np.bitwise_or(UNSA, DNU)

    out_model = undersampling_correction(ramp_model, signal_threshold)
    out_gdq = out_model.groupdq

    npt.assert_array_equal(out_gdq, true_out_gdq)


def test_pix_2():
    """
    Tests groups having data exceeding the signal threshold, Some groups are
    already flagged as DO_NOT_USE; they will not be checked for UC,  Other
    groups will have 'DNU'+'UNSA' added to their GROUPDQ, as will all later
    groups.
    """
    ngroups, nints, nrows, ncols = set_scalars()
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    signal_threshold = 20000.

    # Populate SCI and GROUPDQ arrays.
    ramp_model.data[0, 1, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.data[0, 2, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)
    ramp_model.data[0, 3, 0, 0] = np.array((signal_threshold + 100.), dtype=np.float32)

    ramp_model.groupdq[0, 1, 0, 0] = DNU  # should not get UNSA
    ramp_model.groupdq[0, 2, 0, 0] = np.bitwise_or(ADFL, DNU)  # should not get UNSA
    ramp_model.groupdq[0, 3, 0, 0] = ADFL  # should get UNSA + DNU

    true_out_gdq = ramp_model.groupdq.copy()
    true_out_gdq[0, 3, 0, 0] = np.bitwise_or(np.bitwise_or(DNU, UNSA), ADFL)
    true_out_gdq[0, 4:, 0, 0] = np.bitwise_or(DNU, UNSA)

    out_model = undersampling_correction(ramp_model, signal_threshold)
    out_gdq = out_model.groupdq

    npt.assert_array_equal(out_gdq, true_out_gdq)


def test_too_few_groups():
    """
    Test that processing for datasets having too few (<3) groups per integration
    are skipped.
    """
    ngroups, nints, nrows, ncols = set_scalars()
    ngroups = 2

    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    ramp_model.data[0, :, 0, 0] = 20000.
    sig_thresh = 100.

    result = UndersamplingCorrectionStep.call(ramp_model, skip=False,
                                              signal_threshold=sig_thresh)
    status = result.meta.cal_step.undersampling_correction

    npt.assert_string_equal(status, "SKIPPED")


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
    ERR extensions, and create datamodels for the ramp. All groups
    in all arrays have initial values set to 0.
    """
    err = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint8)

    # Create and populate ramp model
    ramp_model = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq)
    ramp_model.meta.instrument.name = 'NIRISS'
    ramp_model.meta.instrument.detector = 'NIS'

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
