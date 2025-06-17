import numpy as np

from stdatamodels.jwst.datamodels import RampModel, dqflags

from jwst.rscd.rscd_sub import correction_skip_groups


def test_rscd_baseline_set_groupdq():
    """
    For a 2 integration exposure, test if the rscd code sets the groupdq flag on
    the n groups to 'do_not_use' for the 2nd integration and did not change the
    groupdq flags in the 1st integration
    """

    exposure = {"integration_start": None, "integration_end": None, "ngroups": 10, "nints": 2}
    # size of integration
    ngroups = exposure["ngroups"]
    nints = exposure["nints"]

    xsize = 10
    ysize = 10

    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0, dtype=np.float32)
    groupdq = np.zeros(csize, dtype=np.uint8)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)
    dm_ramp.meta.exposure._instance.update(exposure)

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


def test_rscd_baseline_too_few_groups():
    """
    Test that the baseline algorithm is skipped if too few groups are present
    in the integrations.
    """

    # size of exposure
    xsize = 10
    ysize = 10

    exposure = {"integration_start": None, "integration_end": None, "ngroups": 3, "nints": 2}
    # size of integration
    ngroups = exposure["ngroups"]
    nints = exposure["nints"]

    xsize = 10
    ysize = 10

    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0, dtype=np.float32)
    groupdq = np.zeros(csize, dtype=np.uint8)

    # create a JWST datamodel for MIRI data on a copy (the copy is created at the step script)
    dm_ramp = RampModel(data=data, groupdq=groupdq)
    dm_ramp.meta.exposure._instance.update(exposure)

    # get the number of groups to flag
    nflag = 3

    # run the RSCD baseline correction step
    dm_ramp_rscd = correction_skip_groups(dm_ramp.copy(), nflag)

    # test that the groupdq flags are not changed for any integration
    dq_diff = dm_ramp_rscd.groupdq[:, :, :, :] - dm_ramp.groupdq[:, :, :, :]
    np.testing.assert_array_equal(
        np.full((nints, ngroups, ysize, xsize), 0, dtype=int),
        dq_diff,
        err_msg="groupdq flags changed when " + "not enough groups are present",
    )


def test_rscd_tso():
    """
    The RSCD correction is generally skipped for TSO data, but some users
    have been running it on TSO data. So this test was added.
    Test for TSO segmented data that the correct groups are flagged as 'DO_NOT_USE'
    for integration 2 and higher. Set up the segmented data so the segment
    is for integration 25 to 50. A rscd correction should be applied to all
    the data.

    """
    exposure = {"integration_start": 25, "integration_end": 50, "ngroups": 8, "nints": 50}

    xsize = 10
    ysize = 10
    ngroups = exposure["ngroups"]
    seg_nints = exposure["integration_end"] - exposure["integration_start"] + 1
    input_model = RampModel((seg_nints, exposure["ngroups"], ysize, xsize))

    input_model.groupdq[:, :, :, :] = 0  # initize to 0 - all good
    input_model.meta.exposure._instance.update(exposure)

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
