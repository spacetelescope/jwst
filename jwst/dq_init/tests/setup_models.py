from jwst.datamodels import MIRIRampModel
from jwst.datamodels import ImageModel
from jwst.datamodels import MaskModel
import numpy as np


def make_rawramp(nints,ngroups,ysize,xsize):
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    # pixeldq = np.zeros((ysize,xsize), dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data)
    if (ysize == 1024) and (xsize == 1032):
        dm_ramp.meta.subarray.xstart = 1
        dm_ramp.meta.subarray.xsize = xsize
        dm_ramp.meta.subarray.ystart = 1
        dm_ramp.meta.subarray.ysize = ysize
    else:
        # how to set up for subarray... assume
        # specific subarray values or allow
        # user to specify?
        dm_ramp.meta.subarray.xstart = 1
        dm_ramp.meta.subarray.xsize = xsize
        dm_ramp.meta.subarray.ystart = 1 # 467 for MASK1550
        dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp


def make_rampmodel(nints, ngroups, ysize, xsize):
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

    dm_ramp.meta.instrument.name = 'MIRI'
    dm_ramp.meta.observation.date = '2018-01-01'
    dm_ramp.meta.observation.time = '00:00:00'
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp


def make_imagemodel(ysize, xsize):
    # create the data array
    csize = (ysize, xsize)
    data = np.full(csize, 1.0)
    err = np.zeros(csize)
    dq = np.zeros(csize, dtype=int)
    im = ImageModel(data=data, err=err, dq=dq)

    return im


def make_maskmodel(ysize, xsize):
    # create a mask model for the dq_init step
    csize = (ysize, xsize)
    dq = np.zeros(csize, dtype=int)
    # how do we define a dq_def extension?
    mask = MaskModel()

    dqdef = [(0, 1, 'DO_NOT_USE', 'Bad Pixel do not use'),
             (1, 2, 'DEAD', 'Dead Pixel'),
             (2, 4, 'HOT', 'Hot pixel'),
             (3, 8, 'UNRELIABLE_SLOPE', 'Large slope variance'),
             (4, 16, 'RC', 'RC pixel'),
             (5, 32, 'REFERENCE_PIXEL', 'Reference Pixel')]

    dq_def = np.array((dqdef), dtype=mask.dq_def.dtype)
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = 'MIRI'
    ref_data.meta.subarray.xstart = 1
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = 1
    ref_data.meta.subarray.ysize = ysize

    return ref_data
