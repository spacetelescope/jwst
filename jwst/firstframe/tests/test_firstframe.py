import numpy as np
import pytest

from jwst.datamodels.miri_ramp import MIRIRampModel
from jwst.datamodels import dqflags
from jwst.firstframe.firstframe_sub import do_correction

def test_firstframe_set_groupdq():
    """ 
    Test if the firstframe code set the groupdq flag on the first
    group to 'do_not_use'
    """

    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, xsize, ysize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data, groupdq=groupdq)

    # run the first frame correction step
    dm_ramp_firstframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag
    dq_diff = dm_ramp_firstframe.groupdq[0,0,:,:] - dm_ramp.groupdq[0,0,:,:]
    
    np.testing.assert_array_equal(np.full((xsize,ysize),
                                          dqflags.group['DO_NOT_USE'],
                                          dtype=int),
                                  dq_diff,
                                  err_msg='Diff in groupdq flags is not ' \
                                           + 'equal to the DO_NOT_USE flag')

    # test that the groupdq flags are not changed for the rest of the groups
    dq_diff = (dm_ramp_firstframe.groupdq[0,1:ngroups,:,:]
               - dm_ramp.groupdq[0,1:ngroups,:,:])
    np.testing.assert_array_equal(np.full((ngroups-1,xsize,ysize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='n >= 2 groupdq flags changes ' \
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
    csize = (1, ngroups, xsize, ysize)
    data = np.full(csize, 1.0)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data, groupdq=groupdq)

    # run the first frame correction step
    dm_ramp_firstframe = do_correction(dm_ramp)

    # check that the difference in the groupdq flags is equal to
    #   the 'do_not_use' flag

    dq_diff = dm_ramp_firstframe.groupdq[0,0,:,:] - dm_ramp.groupdq[0,0,:,:]
    
    np.testing.assert_array_equal(np.full((xsize,ysize),
                                          0,
                                          dtype=int),
                                  dq_diff,
                                  err_msg='groupdq changed for single group ' \
                                           + 'when it should not')
    
if __name__ == '__main__':
    test_firstframe_set_groupdq()
    test_firstframe_single_group()
