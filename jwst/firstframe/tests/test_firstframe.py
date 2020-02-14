import numpy as np

from jwst.datamodels import RampModel
from jwst.datamodels import dqflags
from jwst.firstframe.firstframe_sub import do_correction
from jwst.firstframe import FirstFrameStep


def test_firstframe_set_groupdq():
    """
    Test if the firstframe code set the groupdq flag on the first
    group to 'do_not_use' for 5 integrations
    """

    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # run the first frame correction step
    dm_ramp_firstframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag
    dq_diff = dm_ramp_firstframe.groupdq[0, 0, :, :] - dm_ramp.groupdq[0, 0, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          dqflags.group['DO_NOT_USE'],
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                  + 'equal to the DO_NOT_USE flag')

    # test that the groupdq flags are not changed for the rest of the groups
    dq_diff = (dm_ramp_firstframe.groupdq[0, 1:ngroups, :, :]
               - dm_ramp.groupdq[0, 1:ngroups, :, :])
    np.testing.assert_array_equal(np.full((ngroups-1, ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='n >= 2 groupdq flags changes '
                                  + 'and they should not be')


def test_firstframe_single_group():
    """
    Test that the firstframe code does nothing when passed a single
    group integration
    """

    # size of integration
    ngroups = 1
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # run the first frame correction step
    dm_ramp_firstframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag

    dq_diff = dm_ramp_firstframe.groupdq[0, 0, :, :] - dm_ramp.groupdq[0, 0, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq changed for single group '
                                  + 'when it should not')


def test_firstframe_add1_groupdq():
    """
    Test if the firstframe code set the groupdq flag on the first
    group to 'do_not_use' by adding 1 to the flag, not overwritting to 1
    """

    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # set a flag in the groupdq, first frame
    dm_ramp.groupdq[0, 0, 500:510, 500:510] = 4

    # run the first frame correction step
    dm_ramp_firstframe = do_correction(dm_ramp)

    # test if pixels in groupdq were incremented in value by 1
    assert(dm_ramp_firstframe.groupdq[0, 0, 505, 505] == 5)


def test_firstframe_3groups():
    """
    Test if the firstframe code set the groupdq flag on the first
    group to 'do_not_use' or left it as is, which it should do for 3 frames
    """

    # size of integration
    ngroups = 3
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # run the first frame correction step
    dm_ramp_firstframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to

    #  0 since the step should not run if frames 3 or fewer
    dq_diff = dm_ramp_firstframe.groupdq[0, 0, :, :] - dm_ramp.groupdq[0, 0, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                  + 'equal to 0')


def test_nircam():
    # test that the code skips processing for a NIR instrument

    # size of integration
    ngroups = 3
    xsize = 2048
    ysize = 2048

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    dm_ramp.meta.instrument.name = 'NIRCAM'
    dm_ramp.meta.instrument.detector = 'NRCA1'

    # run the first frame correction step
    dm_ramp_firstframe = FirstFrameStep.call(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   the 0 since the step should not run for NIR data
    dq_diff = dm_ramp_firstframe.groupdq[0, 0, :, :] - dm_ramp.groupdq[0, 0, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                          + 'equal to 0')


def test_miri():
    # test that the code runs given MIRI as an instrument

    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    dm_ramp.meta.instrument.name = 'MIRI'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'

    # run the first frame correction step
    dm_ramp_firstframe = FirstFrameStep.call(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #  DO_NOT_USE flag
    dq_diff = dm_ramp_firstframe.groupdq[0, 0, :, :] - dm_ramp.groupdq[0, 0, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          dqflags.group['DO_NOT_USE'],
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                          + 'equal to DO_NOT_USE')
