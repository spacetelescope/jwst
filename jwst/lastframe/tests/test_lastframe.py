import numpy as np

from jwst.datamodels import RampModel
from jwst.datamodels import dqflags
from jwst.lastframe.lastframe_sub import do_correction
from jwst.lastframe import LastFrameStep


def test_lastframe_set_groupdq():
    """
    Test if the lastframe code set the groupdq flag on the last
    group to 'do_not_use'. For ngroups >= 3, lastframe should be flagged.
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

    # run the last frame correction step
    dm_ramp_lastframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag
    dq_diff = dm_ramp_lastframe.groupdq[0, ngroups - 1, :, :] - dm_ramp.groupdq[0, ngroups - 1, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          dqflags.group['DO_NOT_USE'],
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                  + 'equal to the DO_NOT_USE flag')

    # test that the groupdq flags are not changed for the rest of the groups
    dq_diff = (dm_ramp_lastframe.groupdq[0, 0:ngroups - 2, :, :]
               - dm_ramp.groupdq[0, 0:ngroups - 2, :, :])
    np.testing.assert_array_equal(np.full((ngroups - 2, ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='n <= ngroups-2 groupdq flags changes '
                                  + 'and they should not be')


def test_lastframe_ngroup2():
    """
    Test if the lastframe code set the groupdq flag on the last
    group to 'do_not_use'. Lastframe should not be flagged with ngroups=2.
    """

    # size of integration
    ngroups = 2
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, groupdq=groupdq)

    # run the last frame correction step
    dm_ramp_lastframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #  zero
    dq_diff = dm_ramp_lastframe.groupdq[0, ngroups - 1, :, :] - dm_ramp.groupdq[0, ngroups - 1, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq flag changed '
                                          + 'when it should not')


def test_lastframe_single_group():
    """
    Test that the lastframe code does nothing when passed a single
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

    # run the last frame correction step
    dm_ramp_lastframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    # zero

    dq_diff = dm_ramp_lastframe.groupdq[0, ngroups - 1, :, :] - dm_ramp.groupdq[0, ngroups - 1, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq changed for single group '
                                  + 'when it should not')


def test_lastframe_add1_groupdq():
    """
    Test if the lastframe code set the groupdq flag on the first
    group to 'do_not_use' by adding 1 to the flag, not overwriting to 1
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

    # set a flag in the groupdq, last frame
    dm_ramp.groupdq[0, ngroups - 1, 500:510, 500:510] = 4

    # run the last frame correction step
    dm_ramp_lastframe = do_correction(dm_ramp)

    # test if pixels in groupdq were incremented in value by 1
    assert (dm_ramp_lastframe.groupdq[0, ngroups - 1, 505, 505] == 5)


def test_nircam():
    # test that the code skips processing for a NIR instrument

    # size of integration
    ngroups = 5
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

    # run the last frame correction step
    dm_ramp_lastframe = LastFrameStep.call(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   0 since the step should not run for NIR data
    dq_diff = dm_ramp_lastframe.groupdq[0, ngroups - 1, :, :] - dm_ramp.groupdq[0, ngroups - 1, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                          + 'equal to 0')


def test_miri():
    # test that the code chooses to process a MIRI file

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

    # run the last frame correction step
    dm_ramp_lastframe = LastFrameStep.call(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   DO_NOT_USE flag
    dq_diff = dm_ramp_lastframe.groupdq[0, ngroups - 1, :, :] - dm_ramp.groupdq[0, ngroups - 1, :, :]

    np.testing.assert_array_equal(np.full((ysize, xsize),
                                          dqflags.group['DO_NOT_USE'],
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not '
                                          + 'equal to DO_NOT_USE')
