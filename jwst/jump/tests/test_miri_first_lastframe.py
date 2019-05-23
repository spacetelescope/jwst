"""
Unit test for jump step to test not using MIRI first frame or last frame if set to DO_NOT_USE
"""
import pytest
import numpy as np
from jwst.jump import JumpStep
from jwst import datamodels
from jwst.datamodels import dqflags, MIRIRampModel
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


@pytest.fixture(scope='function')
def science_image():
    """Generate science image"""
    # size of integration
    nints = 1
    ngroups = 7
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    # create a JWST datamodel for MIRI data
    im = MIRIRampModel((nints, ngroups, ysize, xsize))
    im.meta.instrument.detector = 'MIRIMAGE'
    im.meta.instrument.name = 'MIRI'
    im.meta.subarray.name = 'FULL'
    im.meta.subarray.xstart = 1
    im.meta.subarray.ystart = 1
    im.meta.subarray.xsize = 1032
    im.meta.subarray.ysize = 1024

    im.meta.observation.date = '2019-02-27'
    im.meta.observation.time = '13:37:18.548'
    im.meta.date = '2019-02-27T13:37:18.548'
    im.meta.exposure.readpatt = 'FAST'
    im.meta.exposure.nframes = ngroups
    im.meta.exposure.group_time = 2.77504

    # set up the data so that if the first and last group are used in jump
    # detection it would cause a jump to be detected between group 1-2
    # and group 6-7
    # set group 1 to be 10,000
    im.data[0, 0,:,:] = 10000.0
    # set groups 1,2 - to be around 30,000
    im.data[0, 1,:,:] = 30000.0
    im.data[0, 2,:,:] = 30020.0
    # set up a jump
    im.data[0, 3,:,:] = 40040.0
    im.data[0, 4,:,:] = 40060.0
    im.data[0, 5,:,:] = 40080.0

    # set group 6 to be 50,000
    im.data[0, 6,:,:] = 50000.0

    # set the first and last group to DO_NOTUSE
    im.groupdq[0, 0,:,:] = dqflags.group['DO_NOT_USE']
    im.groupdq[0, -1,:,:] = dqflags.group['DO_NOT_USE']

    return im


def test_miri_first_lastframe(_jail, science_image):

    im = datamodels.open(science_image)
    xsize = im.meta.subarray.xsize
    ysize = im.meta.subarray.ysize

    collect_pipeline_cfgs('./config')
    result = JumpStep.call(im, config_file='config/jump.cfg')

    first_group = result.groupdq[0,0,:,:]
    last_group = result.groupdq[0,6,:,:]
    jump_group = result.groupdq[0,3,:,:]

    do_not_use_array = np.zeros((ysize,xsize))
    jump_array = np.zeros((ysize,xsize))
    do_not_use_array[:,:] = dqflags.group['DO_NOT_USE']
    jump_array[:,:] = dqflags.group['JUMP_DET']

    assert np.array_equal(first_group, do_not_use_array)
    assert np.array_equal(last_group, do_not_use_array)
    assert np.array_equal(jump_group, jump_array)
