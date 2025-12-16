import numpy as np
from stdatamodels.jwst.datamodels import dqflags

from jwst.rscd.rscd_sub import correction_skip_groups


def test_rscd_baseline_set_groupdq(create_miri_model):
    """
    For a 2 integration exposure, test if the rscd code sets the correct
    groupdq flag on  the n groups to 'do_not_use' for each integration
    """
    ngroups = 10
    nints = 2
    xsize = 10
    ysize = 10
    dm_ramp = create_miri_model(nints, ngroups, ysize, xsize)

    # set the number of groups to flag in int 1 and int 2
    nflag_int1 = 1
    nflag_int2 = 3

    # run the RSCD baseline correction step on a copy (the copy is created at the step script)
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag_int1, nflag_int2)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag for the 2nd integration
    dq_diff = dm_ramp_rscd.groupdq[1, 0:nflag_int2, :, :] - dm_ramp.groupdq[1, 0:nflag_int2, :, :]
    np.testing.assert_array_equal(
        np.full((nflag_int2, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 2nd integration
    dq_diff = (
        dm_ramp_rscd.groupdq[1, nflag_int2:ngroups, :, :]
        - dm_ramp.groupdq[1, nflag_int2:ngroups, :, :]
    )
    np.testing.assert_array_equal(
        np.full((ngroups - nflag_int2, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )

    # test that the groupdq flags for first integration are correct
    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag for the 1st integration
    dq_diff = dm_ramp_rscd.groupdq[0, 0:nflag_int1, :, :] - dm_ramp.groupdq[0, 0:nflag_int1, :, :]
    np.testing.assert_array_equal(
        np.full((nflag_int1, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 1st integration
    dq_diff = (
        dm_ramp_rscd.groupdq[0, nflag_int1:ngroups, :, :]
        - dm_ramp.groupdq[0, nflag_int1:ngroups, :, :]
    )
    np.testing.assert_array_equal(
        np.full((ngroups - nflag_int1, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )


def test_rscd_baseline_set_groupdq_test2(create_miri_model):
    """
    For a 2 integration exposure, set the number of groups to reject
    more than the would leave us with 3 good groups. Test that the number
    actually rejected is what will give us 3 remaining groups
    """
    ngroups = 10
    nints = 2
    xsize = 10
    ysize = 10
    dm_ramp = create_miri_model(nints, ngroups, ysize, xsize)

    # set the number of groups to flag in int 1 and int 2
    nflag_int1 = 8
    nflag_int2 = 10

    # run the RSCD baseline correction step on a copy (the copy is created at the step script)
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag_int1, nflag_int2)

    # what should happen
    nflag_int2_used = 7
    nflag_int1_used = 7
    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag for the 2nd integration
    dq_diff = (
        dm_ramp_rscd.groupdq[1, 0:nflag_int2_used, :, :]
        - dm_ramp.groupdq[1, 0:nflag_int2_used, :, :]
    )
    np.testing.assert_array_equal(
        np.full((nflag_int2_used, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 2nd integration
    dq_diff = (
        dm_ramp_rscd.groupdq[1, nflag_int2_used:ngroups, :, :]
        - dm_ramp.groupdq[1, nflag_int2_used:ngroups, :, :]
    )
    np.testing.assert_array_equal(
        np.full((ngroups - nflag_int2_used, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )

    # test that the groupdq flags for first integration are correct
    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag for the 1st integration
    dq_diff = (
        dm_ramp_rscd.groupdq[0, 0:nflag_int1_used, :, :]
        - dm_ramp.groupdq[0, 0:nflag_int1_used, :, :]
    )
    np.testing.assert_array_equal(
        np.full((nflag_int1_used, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 1st integration
    dq_diff = (
        dm_ramp_rscd.groupdq[0, nflag_int1_used:ngroups, :, :]
        - dm_ramp.groupdq[0, nflag_int1_used:ngroups, :, :]
    )
    np.testing.assert_array_equal(
        np.full((ngroups - nflag_int1_used, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )


def test_rscd_baseline_too_few_groups(create_miri_model):
    """
    Test that the baseline algorithm is skipped if too few groups are present
    in the integrations.
    """
    dm_ramp = create_miri_model(ngroups=2)

    # get the number of groups to flag
    nflag_int1 = 3
    nflag_int2 = 3
    # run the RSCD baseline correction step
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag_int1, nflag_int2)

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

    # get the number of groups to flag in int 1 and int 2
    nflag_int1 = 1
    nflag_int2 = 4

    # run the RSCD baseline correction step on a copy (the copy is created at the step script)
    ramp_rscd = correction_skip_groups(input_model.copy(), nflag_int1, nflag_int2)

    # check that the difference in the groupdq flags is equal to
    #  the 'DO_NOT_USE' flag for the 1st integration in the segment
    # which is actually the 25th integration
    dq_diff = ramp_rscd.groupdq[0, 0:nflag_int2, :, :] - input_model.groupdq[0, 0:nflag_int2, :, :]

    np.testing.assert_array_equal(
        np.full((nflag_int2, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int),
        dq_diff,
        err_msg="Diff in groupdq flags is not " + "equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 1st  integration in the segment
    dq_diff = (
        ramp_rscd.groupdq[0, nflag_int2:ngroups, :, :]
        - input_model.groupdq[0, nflag_int2:ngroups, :, :]
    )
    np.testing.assert_array_equal(
        np.full((ngroups - nflag_int2, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed after " + "maximum requested",
    )


def test_rscd_addto_groupdq(create_miri_model):
    """
    Test that the groupdq flag for flagged frames as 'DO_NOT_USE' has
    1 added to the flag and not overwriting to 1.
    """
    ngroups = 10
    nints = 1
    xsize = 10
    ysize = 10
    dm_ramp = create_miri_model(nints, ngroups, ysize, xsize)
    dm_ramp.groupdq[:, 0, 5, 5] = 4
    # set the number of groups to flag in int 1 and int 2
    nflag_int1 = 1
    nflag_int2 = 1

    # run the RSCD baseline correction step on a copy (the copy is created at the step script
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag_int1, nflag_int2)

    # test if pixels in groupdq were incremented in value by 1
    assert dm_ramp_rscd.groupdq[0, 0, 5, 5] == 5


def test_rscd_bright_use_group1(create_miri_model):
    """
    Test the rscd code when bright_use_group1 is set to True.

    The groupdq flag for group 1 should not be set to DO_NOT_USE for the pixels that saturate
    the 3rd group.  Otherwise, all other pixels should have their group1 groupdq
    flags set to DO_NOT_USE.
    """

    # size of ramp
    nints = 2
    ngroups = 5
    xsize = 40
    ysize = 40

    dm_ramp = create_miri_model(nints, ngroups, ysize, xsize)
    # create the data and groupdq arrays
    dm_ramp.data[:, :, :, :] = 1

    # set a fraction of the pixels to saturate in between the 2nd and 3rd groups
    dm_ramp.groupdq[:, 2, 0:10, :] = dqflags.group["SATURATED"]

    # set a fraction of the pixels to saturate in between the 1st and 2nd groups
    dm_ramp.groupdq[:, 2, 20:30, :] = dqflags.group["SATURATED"]  # saturates on group 3
    dm_ramp.groupdq[:, 1, 20:30, :] = dqflags.group["SATURATED"]  # saturates on group 2
    #  Set a pixels to also be saturated on the first group. These will be flagged
    # by rscd correction.

    dm_ramp.groupdq[:, 0, 25:30, :] = dqflags.group["SATURATED"]  # saturates on group 1

    dm_ramp.pixeldq[:, :] = 0  # initialize to zero

    # check for a log message
    n_kept = 15 * xsize * nints
    n_kept_int = 15 * xsize

    nflag_int1 = 1
    nflag_int2 = 1
    dm_ramp_rscd = correction_skip_groups(
        dm_ramp.copy(), nflag_int1, nflag_int2, bright_use_group1=True
    )
    number_kept1 = dm_ramp_rscd.meta.rscd.keep_bright_firstgroup_int1
    assert number_kept1 == n_kept_int

    number_kept2 = dm_ramp_rscd.meta.rscd.keep_bright_firstgroup_int2p
    assert number_kept2 == n_kept_int

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag
    dq_diff = dm_ramp_rscd.groupdq[:, 0, :, :] - dm_ramp.groupdq[:, 0, :, :]

    expected_diff = np.full((nints, ysize, xsize), dqflags.group["DO_NOT_USE"], dtype=int)
    expected_diff[:, 0:10, :] = 0
    expected_diff[:, 20:25, :] = 0

    np.testing.assert_array_equal(
        expected_diff,
        dq_diff,
        err_msg="Diff in groupdq flags is not equal to the DO_NOT_USE flag",
    )

    # test that the groupdq flags are not changed for the rest of the groups
    dq_diff = dm_ramp_rscd.groupdq[1, 1:ngroups, :, :] - dm_ramp.groupdq[1, 1:ngroups, :, :]
    index = np.where(dq_diff != 0)

    np.testing.assert_array_equal(
        np.full((ngroups - 1, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="n >= 2 groupdq flags changes and they should not be",
    )

    # also check that pixel dq flags are set for the pixels that were kept
    dq_diff = dm_ramp_rscd.pixeldq - dm_ramp.pixeldq
    expected_diff = np.full((ysize, xsize), 0, dtype=int)
    expected_diff[0:10, :] = dqflags.pixel["FLUX_ESTIMATED"]
    expected_diff[20:25, :] = dqflags.pixel["FLUX_ESTIMATED"]
    np.testing.assert_array_equal(
        expected_diff,
        dq_diff,
        err_msg="pixeldq flags do not have expected values",
    )
