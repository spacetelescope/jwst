import numpy as np

from stdatamodels.jwst.datamodels import RampModel, dqflags

from jwst.rscd.rscd_sub import correction_skip_groups


def test_rscd_baseline_set_groupdq():
    """
    For a 2 integration exposure, test if the rscd code sets the groupdq flag on
    the n groups to 'do_not_use' for the 2nd integration and did not change the
    groupdq flags in the 1st integration
    """

    # size of integration
    ngroups = 10
    xsize = 10
    ysize = 10

    # create the data and groupdq arrays
    csize = (2, ngroups, ysize, xsize)
    data = np.full(csize, 1.0, dtype=np.float32)
    groupdq = np.zeros(csize, dtype=np.uint8)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # get the number of groups to flag
    nflag = 3

    # run the RSCD baseline correction step
    dm_ramp_rscd = correction_skip_groups(dm_ramp, nflag)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag for the 2nd integration
    dq_diff = (dm_ramp_rscd.groupdq[1, 0:nflag, :, :]
               - dm_ramp.groupdq[1, 0:nflag, :, :])
    np.testing.assert_array_equal(np.full((nflag, ysize, xsize),
                                          dqflags.group['DO_NOT_USE'],
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                  + 'equal to the DO_NOT_USE flag')

    # test that the groupdq flags are not changed for the rest of the groups
    # in the 2nd integration
    dq_diff = (dm_ramp_rscd.groupdq[1, nflag:ngroups, :, :]
               - dm_ramp.groupdq[1, nflag:ngroups, :, :])
    np.testing.assert_array_equal(np.full((ngroups - nflag, ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq flags changed after '
                                  + 'maximum requested')

    # test that the groupdq flags are not changed for the 1st integration
    dq_diff = (dm_ramp_rscd.groupdq[0, :, :, :]
               - dm_ramp.groupdq[0, :, :, :])
    np.testing.assert_array_equal(np.full((ngroups, ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq flags changed in 1st '
                                  + 'integration this should not happen')


def test_rscd_baseline_too_few_groups():
    """
    Test that the baseline algorithm is skipped if too few groups are present
    in the integrations.
    """

    # size of exposure
    nints = 2
    ngroups = 3
    xsize = 10
    ysize = 10

    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0, dtype=np.float32)
    groupdq = np.zeros(csize, dtype=np.uint8)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # get the number of groups to flag
    nflag = 3

    # run the RSCD baseline correction step
    dm_ramp_rscd = correction_skip_groups(dm_ramp, nflag)

    # test that the groupdq flags are not changed for any integration
    dq_diff = (dm_ramp_rscd.groupdq[:, :, :, :]
               - dm_ramp.groupdq[:, :, :, :])
    np.testing.assert_array_equal(np.full((nints, ngroups, ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq flags changed when '
                                  + 'not enough groups are present')
