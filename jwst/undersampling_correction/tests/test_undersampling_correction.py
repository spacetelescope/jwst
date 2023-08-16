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
DROU = test_dq_flags["DROPOUT"]
UNSA_DNU = UNSA + DNU


def test_pix_0():
    """
    Having all data in ramp below the signal threshold, the only non-GOOD
    groups in the output GROUPDQ should be those DNU propagated from the input.
    """
    ngroups, nints, nrows, ncols = 10, 1, 1, 1
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
    Tests groups whose data exceeds the signal threshold; 1 group is already
    flagged as DNU from a previous calibration step, and 1 group is GOOD.
    Also tests groups whose data does not exceed the signal threshold;
    similarly 1 group is already flagged as DNU from a previous calibration
    step, and 1 is GOOD. All data beyond the first exceedance are also
    flagged as UNDERSAMP and DNU.
    """
    ngroups, nints, nrows, ncols = 10, 1, 1, 1
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    signal_threshold = 20000.

    # Populate SCI and GROUPDQ arrays.
    ramp_model.data[0, 1, 0, 0] = np.array((0.5 * signal_threshold), dtype=np.float32)

    ramp_model.data[0, 2, 0, 0] = np.array((0.8 * signal_threshold), dtype=np.float32)
    ramp_model.groupdq[0, 2, 0, 0] = DNU  # should not get UNSA, not an exceedance

    ramp_model.data[0, 3, 0, 0] = np.array((signal_threshold + 5000.), dtype=np.float32)
    ramp_model.groupdq[0, 3, 0, 0] = DNU  # should not get UNSA, although exceedance

    ramp_model.data[0, 4:, 0, 0] = np.array((signal_threshold + 6000.), dtype=np.float32)
    ramp_model.groupdq[0, 4:, 0, 0] = GOOD

    true_out_gdq = ramp_model.groupdq.copy()
    true_out_gdq[0, 2, 0, 0] = DNU
    true_out_gdq[0, 3, 0, 0] = DNU
    true_out_gdq[0, 4:, 0, 0] = UNSA_DNU

    out_model = undersampling_correction(ramp_model, signal_threshold)
    out_gdq = out_model.groupdq

    npt.assert_array_equal(out_gdq, true_out_gdq)


def test_too_few_groups():
    """
    Test that processing for datasets having too few (<3) groups per integration
    are skipped.
    """
    ngroups, nints, nrows, ncols = 2, 1, 1, 1
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    ramp_model.data[0, :, 0, 0] = 20000.
    sig_thresh = 100.

    result = UndersamplingCorrectionStep.call(ramp_model, skip=False,
                                              signal_threshold=sig_thresh)
    status = result.meta.cal_step.undersampling_correction

    npt.assert_string_equal(status, "SKIPPED")


def test_flag_neighbors():
    """
    Test flagging of 4 nearest neighbors of exceedances. Tests pixels on
    array edges, Tests exclusion of groups previously flagged as DO_NOT_USE.
    """
    ngroups, nints, nrows, ncols = 6, 1, 4, 3
    ramp_model, pixdq, groupdq, err = create_mod_arrays(
        ngroups, nints, nrows, ncols)

    signal_threshold = 4400.

    # Populate pixel-specific SCI and GROUPDQ arrays
    ramp_model.data[0, :, :, :] = \
        np.array([[
            [1900., 2666., 2100.],
            [3865., 2300., 3177.],
            [3832., 3044., 3588.],
            [3799., 3233., 3000.]],

           [[2100., 2866., 2300.],
            [4065., 2500., 3377.],
            [4032., 3244., 3788.],
            [3999., 3433., 3200.]],

           [[2300., 3066., 2500.],
            [4265., 2700., 3577.],
            [4232., 3444., 3988.],
            [4199., 3633., 3400.]],

           [[2500., 3266., 2700.],
            [4465., 2900., 3777.],
            [4432., 3644., 4188.],
            [4399., 3833., 3600.]],

           [[2700., 3466., 2900.],
            [4665., 3100., 3977.],
            [4632., 3844., 4388.],
            [4599., 4033., 3800.]],

           [[2900., 3666., 3100.],
            [4865., 3300., 4177.],
            [4832., 4044., 4588.],
            [4799., 4233., 4000.]]], dtype=np.float32)

    # These group DQ values should propagate unchanged to the output
    ramp_model.groupdq[:, 4, 2, 0] = [DNU]
    ramp_model.groupdq[:, 1, 2, 2] = [DNU]
    ramp_model.groupdq[:, 2, 1, 1] = [DROU + DNU]

    out_model = undersampling_correction(ramp_model, signal_threshold)
    out_gdq = out_model.groupdq

    true_out_gdq = ramp_model.groupdq.copy()
    true_out_gdq[0, :, :, :] = \
        np.array([[
            [0,   0,   0],
            [0,   0,   0],
            [0,   0,   0],
            [0,   0,   0]],

           [[0,   0,   0],
            [0,   0,   0],
            [0,   0,   DNU],
            [0,   0,   0]],

           [[0,   0,   0],
            [0,   9,   0],
            [0,   0,   0],
            [0,   0,   0]],

           [[UNSA_DNU,   0,      0],
            [UNSA_DNU, UNSA_DNU, 0],
            [UNSA_DNU, UNSA_DNU, 0],
            [UNSA_DNU,   0,      0]],

           [[UNSA_DNU,   0,      0],
            [UNSA_DNU, UNSA_DNU, 0],
            [DNU,        0,      0],
            [UNSA_DNU, UNSA_DNU, 0]],

           [[UNSA_DNU,   0,       0],
            [UNSA_DNU, UNSA_DNU, UNSA_DNU],
            [UNSA_DNU, UNSA_DNU, UNSA_DNU],
            [UNSA_DNU, UNSA_DNU, UNSA_DNU]]], dtype=np.uint8)

    npt.assert_array_equal(out_gdq, true_out_gdq)


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
