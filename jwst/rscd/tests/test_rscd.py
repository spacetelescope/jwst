import numpy as np
from stdatamodels.jwst.datamodels import dqflags

from jwst.rscd.rscd_sub import correction_skip_groups


def test_rscd_baseline_set_groupdq(create_miri_model):
    """
    For a 2 integration exposure, test if the rscd code sets the groupdq flag on
    the n groups to 'do_not_use' for the 2nd integration and did not change the
    groupdq flags in the 1st integration
    """
    ngroups = 10
    nints = 2
    xsize = 10
    ysize = 10
    dm_ramp = create_miri_model(nints, ngroups, ysize, xsize)

    # get the number of groups to flag
    nflag = 3

    # run the RSCD baseline correction step on a copy (the copy is created at the step script)
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag for the 2nd integration
    dq_diff = dm_ramp_rscd.groupdq[1, 0:nflag, :, :] - dm_ramp.groupdq[1, 0:nflag, :, :]
    np.testing.assert_array_equal(
        np.full((nflag, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 2nd integration
    dq_diff = dm_ramp_rscd.groupdq[1, nflag:ngroups, :, :] - dm_ramp.groupdq[1, nflag:ngroups, :, :]
    np.testing.assert_array_equal(
        np.full((ngroups - nflag, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )

    # test that the groupdq flags are not changed for the 1st integration
    dq_diff = dm_ramp_rscd.groupdq[0, :, :, :] - dm_ramp.groupdq[0, :, :, :]
    np.testing.assert_array_equal(
        np.full((ngroups, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed in 1st " + "integration this should not happen",
    )


def test_rscd_baseline_too_few_groups(create_miri_model):
    """
    Test that the baseline algorithm is skipped if too few groups are present
    in the integrations.
    """
    dm_ramp = create_miri_model(ngroups=3)

    # get the number of groups to flag
    nflag = 3

    # run the RSCD baseline correction step
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag)

    # test that the groupdq flags are not changed for any integration
    dq_diff = dm_ramp_rscd.groupdq[:, :, :, :] - dm_ramp.groupdq[:, :, :, :]
    np.testing.assert_array_equal(
        dq_diff,
        0,
        err_msg="groupdq flags changed when not enough groups are present",
    )


def test_rscd_tso(create_miri_model):
    """
    The RSCD correction is generally skipped for TSO data, but some users
    have been running it on TSO data. So this test was added.
    Test for TSO segmented data that the correct groups are flagged as 'DO_NOT_USE'
    for integration 2 and higher. Set up the segmented data so the segment
    is for integration 25 to 50. A rscd correction should be applied to all
    the data.

    """
    xsize = 10
    ysize = 10
    ngroups = 8
    int_start = 25
    int_end = 50
    seg_nints = int_end - int_start + 1
    input_model = create_miri_model(nints=seg_nints, ngroups=8, ysize=ysize, xsize=xsize)
    input_model.meta.exposure.integration_start = int_start
    input_model.meta.exposure.integration_end = int_end

    # get the number of groups to flag
    nflag = 4

    # run the RSCD baseline correction step on a copy (the copy is created at the step script)
    ramp_rscd = correction_skip_groups(input_model.copy(), nflag)

    # check that the difference in the groupdq flags is equal to
    #  the 'DO_NOT_USE' flag for the 1st integration in the segment
    # which is actually the 25th integration
    dq_diff = ramp_rscd.groupdq[0, 0:nflag, :, :] - input_model.groupdq[0, 0:nflag, :, :]

    np.testing.assert_array_equal(
        np.full((nflag, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 1st  integration in the segment
    dq_diff = (
        ramp_rscd.groupdq[0, nflag:ngroups, :, :] - input_model.groupdq[0, nflag:ngroups, :, :]
    )
    np.testing.assert_array_equal(
        np.full((ngroups - nflag, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )
