import numpy as np

from stdatamodels.jwst.datamodels import RampModel
from stdatamodels.jwst.datamodels import dqflags

from jwst.charge_migration.charge_migration import charge_migration
from jwst.charge_migration.charge_migration_step import ChargeMigrationStep

import numpy.testing as npt


test_dq_flags = dqflags.pixel
GOOD = test_dq_flags["GOOD"]
DNU = test_dq_flags["DO_NOT_USE"]
CHLO = test_dq_flags["CHARGELOSS"]

DROU = test_dq_flags["DROPOUT"]
CHLO_DNU = CHLO + DNU


def test_pix_0():
    """
    Having all data in ramp below the signal threshold, the only non-GOOD
    groups in the output GROUPDQ should be those DNU propagated from the input.
    """
    ngroups, nints, nrows, ncols = 10, 1, 1, 1
    ramp_model, pixdq, groupdq, err = create_mod_arrays(ngroups, nints, nrows, ncols)

    signal_threshold = 30000.0

    # Populate pixel-specific SCI and GROUPDQ arrays
    # Set all SCI to be below the signal threshold, and some input
    # GROUPDQ to be DNU
    ramp_model.data[0, :, 0, 0] = np.array((signal_threshold - 100.0), dtype=np.float32)
    ramp_model.groupdq[0, :, 0, 0] = [GOOD, DNU, GOOD, DNU, DNU, GOOD, GOOD, GOOD, GOOD, DNU]
    in_data = ramp_model.data
    in_gdq = ramp_model.groupdq

    out_model = charge_migration(ramp_model, signal_threshold)
    out_data = out_model.data
    out_gdq = out_model.groupdq

    npt.assert_array_equal(in_data, out_data)
    npt.assert_array_equal(out_gdq, in_gdq)


def test_pix_1():
    """
    Tests groups whose data exceeds the signal threshold; 1 group is already
    flagged as DNU from a previous calibration step, and 1 group is GOOD.
    Also tests groups whose data does not exceed the signal threshold;
    similarly 1 group is already flagged as DNU from a previous calibration
    step, and 1 is GOOD. All data beyond the first exceedance are also
    flagged as CHARGELOSS and DNU.
    """
    ngroups, nints, nrows, ncols = 10, 1, 1, 1
    ramp_model, pixdq, groupdq, err = create_mod_arrays(ngroups, nints, nrows, ncols)

    signal_threshold = 20000.0

    # Populate SCI and GROUPDQ arrays.
    ramp_model.data[0, 1, 0, 0] = np.array((0.5 * signal_threshold), dtype=np.float32)

    ramp_model.data[0, 2, 0, 0] = np.array((0.8 * signal_threshold), dtype=np.float32)
    ramp_model.groupdq[0, 2, 0, 0] = DNU  # should not get CHLO, not an exceedance

    ramp_model.data[0, 3, 0, 0] = np.array((signal_threshold + 5000.0), dtype=np.float32)
    ramp_model.groupdq[0, 3, 0, 0] = DNU  # should not get CHLO, although exceedance

    ramp_model.data[0, 4:, 0, 0] = np.array((signal_threshold + 6000.0), dtype=np.float32)
    ramp_model.groupdq[0, 4:, 0, 0] = GOOD

    true_out_data = ramp_model.data.copy()
    true_out_gdq = ramp_model.groupdq.copy()
    true_out_gdq[0, 2, 0, 0] = DNU
    true_out_gdq[0, 3, 0, 0] = DNU
    true_out_gdq[0, 4:, 0, 0] = CHLO_DNU

    out_model = charge_migration(ramp_model, signal_threshold)

    out_data = out_model.data
    out_gdq = out_model.groupdq

    npt.assert_array_equal(true_out_data, out_data)
    npt.assert_array_equal(out_gdq, true_out_gdq)


def test_pix_2():
    """
    Test a later group being below the threshold.
    """
    ngroups, nints, nrows, ncols = 10, 1, 1, 1
    ramp_model, pixdq, groupdq, err = create_mod_arrays(ngroups, nints, nrows, ncols)

    signal_threshold = 4000.0

    arr = [1000.0, 2000.0, 4005.0, 4500.0, 5000.0, 5500.0, 3500.0, 6000.0, 6500.0, 3700.0]
    ramp_model.data[0, :, 0, 0] = np.array(arr, dtype=np.float32)
    arr = [0, DNU, 0, 0, 0, 0, 0, 0, 0, 0]
    ramp_model.groupdq[0, :, 0, 0] = np.array(arr, dtype=np.uint8)

    out_model = charge_migration(ramp_model, signal_threshold)

    truth_arr = [
        0,
        DNU,
        CHLO_DNU,
        CHLO_DNU,
        CHLO_DNU,
        CHLO_DNU,
        CHLO_DNU,
        CHLO_DNU,
        CHLO_DNU,
        CHLO_DNU,
    ]
    truth_gdq = np.array(truth_arr, dtype=np.uint8)

    npt.assert_array_equal(truth_gdq, out_model.groupdq[0, :, 0, 0])


def nearest_neighbor_base(chg_thresh, pixel):
    """
    Set up ramp array that is 5, 5 with 10 groups.
    The flagging starts in group 3 (zero based) in the pixel tested.
    """
    nints, ngroups, nrows, ncols = 1, 10, 5, 5
    ramp_model, pixdq, groupdq, err = create_mod_arrays(ngroups, nints, nrows, ncols)

    # Set up dummy data
    base = chg_thresh * 0.05
    base_arr = [float(k + 1) * base for k in range(ngroups)]
    for row in range(nrows):
        for col in range(ncols):
            ramp_model.data[0, :, row, col] = np.array(base_arr, dtype=np.float32)

    # Make CHARGELOSS threshold starting at group 3
    in_row, in_col = pixel
    ramp_model.data[0, 3:, in_row, in_col] += chg_thresh

    return ramp_model, pixdq, groupdq, err


def test_nearest_neighbor_1():
    """
    CHARGELOSS center
    The flagging starts in group 3 (zero based) in the pixel tested.
    """
    chg_thresh = 4000.0
    pixel = (2, 2)
    ramp_model, pixdq, groupdq, err = nearest_neighbor_base(chg_thresh, pixel)
    gdq_check = ramp_model.groupdq.copy()
    ngroups = gdq_check.shape[1]

    out_model = charge_migration(ramp_model, chg_thresh)

    check_pattern = [
        [GOOD, GOOD, GOOD, GOOD, GOOD],
        [GOOD, GOOD, CHLO_DNU, GOOD, GOOD],
        [GOOD, CHLO_DNU, CHLO_DNU, CHLO_DNU, GOOD],
        [GOOD, GOOD, CHLO_DNU, GOOD, GOOD],
        [GOOD, GOOD, GOOD, GOOD, GOOD],
    ]
    check = np.array(check_pattern, dtype=gdq_check.dtype)
    for group in range(3, ngroups):
        gdq_check[0, group, :, :] = check

    npt.assert_array_equal(out_model.data, ramp_model.data)
    npt.assert_array_equal(out_model.groupdq, gdq_check)


def test_nearest_neighbor_2():
    """
    CHARGELOSS corner
    The flagging starts in group 3 (zero based) in the pixel tested.
    """
    chg_thresh = 4000.0
    pixel = (0, 0)
    ramp_model, pixdq, groupdq, err = nearest_neighbor_base(chg_thresh, pixel)
    gdq_check = ramp_model.groupdq.copy()
    ngroups = gdq_check.shape[1]

    out_model = charge_migration(ramp_model, chg_thresh)

    check_pattern = [
        [CHLO_DNU, CHLO_DNU, GOOD, GOOD, GOOD],
        [CHLO_DNU, GOOD, GOOD, GOOD, GOOD],
        [GOOD, GOOD, GOOD, GOOD, GOOD],
        [GOOD, GOOD, GOOD, GOOD, GOOD],
        [GOOD, GOOD, GOOD, GOOD, GOOD],
    ]
    check = np.array(check_pattern, dtype=gdq_check.dtype)
    for group in range(3, ngroups):
        gdq_check[0, group, :, :] = check

    npt.assert_array_equal(out_model.data, ramp_model.data)
    npt.assert_array_equal(out_model.groupdq, gdq_check)


def test_nearest_neighbor_3():
    """
    CHARGELOSS Edge
    The flagging starts in group 3 (zero based) in the pixel tested.
    """
    chg_thresh = 4000.0
    pixel = (2, 4)
    ramp_model, pixdq, groupdq, err = nearest_neighbor_base(chg_thresh, pixel)
    gdq_check = ramp_model.groupdq.copy()
    ngroups = gdq_check.shape[1]

    out_model = charge_migration(ramp_model, chg_thresh)

    check_pattern = [
        [GOOD, GOOD, GOOD, GOOD, GOOD],
        [GOOD, GOOD, GOOD, GOOD, CHLO_DNU],
        [GOOD, GOOD, GOOD, CHLO_DNU, CHLO_DNU],
        [GOOD, GOOD, GOOD, GOOD, CHLO_DNU],
        [GOOD, GOOD, GOOD, GOOD, GOOD],
    ]
    check = np.array(check_pattern, dtype=gdq_check.dtype)
    for group in range(3, ngroups):
        gdq_check[0, group, :, :] = check

    npt.assert_array_equal(out_model.data, ramp_model.data)
    npt.assert_array_equal(out_model.groupdq, gdq_check)


def test_too_few_groups():
    """
    Test that processing for datasets having too few (<3) groups per integration
    are skipped.
    """
    ngroups, nints, nrows, ncols = 2, 1, 1, 1
    ramp_model, pixdq, groupdq, err = create_mod_arrays(ngroups, nints, nrows, ncols)

    ramp_model.data[0, :, 0, 0] = 20000.0
    sig_thresh = 100.0

    result = ChargeMigrationStep.call(ramp_model, skip=False, signal_threshold=sig_thresh)
    status = result.meta.cal_step.charge_migration

    npt.assert_string_equal(status, "SKIPPED")


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
    ramp_model.meta.instrument.name = "NIRISS"
    ramp_model.meta.instrument.detector = "NIS"

    ramp_model.meta.subarray.name = "FULL"
    ramp_model.meta.subarray.xstart = 1
    ramp_model.meta.subarray.ystart = 1
    ramp_model.meta.subarray.xsize = ncols
    ramp_model.meta.subarray.ysize = nrows

    ramp_model.meta.exposure.ngroups = ngroups
    ramp_model.meta.exposure.nframes = 1
    ramp_model.meta.exposure.groupgap = 0
    ramp_model.meta.exposure.drop_frames1 = 0
    ramp_model.meta.exposure.type = "NIS_IMAGE"

    return ramp_model, pixdq, gdq, err
