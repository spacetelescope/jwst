import numpy as np
import numpy.testing as npt
from stdatamodels.jwst.datamodels import dqflags

from jwst.charge_migration.charge_migration import flag_charge_migration
from jwst.charge_migration.tests import helpers

GOOD = dqflags.pixel["GOOD"]
DNU = dqflags.pixel["DO_NOT_USE"]
CHLO = dqflags.pixel["CHARGELOSS"]

DROU = dqflags.pixel["DROPOUT"]
CHLO_DNU = CHLO + DNU


def test_pix_0():
    """
    Having all data in ramp below the signal threshold, the only non-GOOD
    groups in the output GROUPDQ should be those DNU propagated from the input.
    """
    shape = (1, 10, 1, 1)
    ramp_model = helpers.niriss_ramp_model(shape)

    signal_threshold = 30000.0

    # Populate pixel-specific SCI and GROUPDQ arrays
    # Set all SCI to be below the signal threshold, and some input
    # GROUPDQ to be DNU
    ramp_model.data[0, :, 0, 0] = np.array((signal_threshold - 100.0), dtype=np.float32)
    ramp_model.groupdq[0, :, 0, 0] = [GOOD, DNU, GOOD, DNU, DNU, GOOD, GOOD, GOOD, GOOD, DNU]
    in_data = ramp_model.data
    in_gdq = ramp_model.groupdq

    out_model = flag_charge_migration(ramp_model, signal_threshold)
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
    shape = (1, 10, 1, 1)
    ramp_model = helpers.niriss_ramp_model(shape)

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

    out_model = flag_charge_migration(ramp_model, signal_threshold)

    out_data = out_model.data
    out_gdq = out_model.groupdq

    npt.assert_array_equal(true_out_data, out_data)
    npt.assert_array_equal(out_gdq, true_out_gdq)


def test_pix_2():
    """
    Test a later group being below the threshold.
    """
    shape = (1, 10, 1, 1)
    ramp_model = helpers.niriss_ramp_model(shape)

    signal_threshold = 4000.0

    arr = [1000.0, 2000.0, 4005.0, 4500.0, 5000.0, 5500.0, 3500.0, 6000.0, 6500.0, 3700.0]
    ramp_model.data[0, :, 0, 0] = np.array(arr, dtype=np.float32)
    arr = [0, DNU, 0, 0, 0, 0, 0, 0, 0, 0]
    ramp_model.groupdq[0, :, 0, 0] = np.array(arr, dtype=np.uint8)

    out_model = flag_charge_migration(ramp_model, signal_threshold)

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


def test_nearest_neighbor_1():
    """
    CHARGELOSS center
    The flagging starts in group 3 (zero based) in the pixel tested.
    """
    chg_thresh = 4000.0
    pixel = (2, 2)
    ramp_model = helpers.niriss_ramp_model_with_signal(chg_thresh, pixel)
    gdq_check = ramp_model.groupdq.copy()
    ngroups = gdq_check.shape[1]

    out_model = flag_charge_migration(ramp_model, chg_thresh)

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
    ramp_model = helpers.niriss_ramp_model_with_signal(chg_thresh, pixel)
    gdq_check = ramp_model.groupdq.copy()
    ngroups = gdq_check.shape[1]

    out_model = flag_charge_migration(ramp_model, chg_thresh)

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
    ramp_model = helpers.niriss_ramp_model_with_signal(chg_thresh, pixel)
    gdq_check = ramp_model.groupdq.copy()
    ngroups = gdq_check.shape[1]

    out_model = flag_charge_migration(ramp_model, chg_thresh)

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
