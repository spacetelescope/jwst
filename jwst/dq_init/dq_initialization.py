#
# Module for handling Data Quality Initialization
#
# This version applies an IRS2-format mask reference file to an IRS2-format
# (3200, 2048) pixeldq array.  Only minor changes were made to the original
# code, in functions is_subarray (use < rather than !=) and check_dimensions
# (pixeldq data type should be uint32).

import logging

import numpy as np

from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def correct_model(input_model, mask_model):
    """DQ Initialize a JWST Model"""

    output_model = do_dqinit (input_model, mask_model)

    return output_model

def do_dqinit (input_model, mask_model):
    """Do the DQ initialization"""

    check_dimensions(input_model)

    output_model = input_model.copy()

    if is_subarray(output_model):
        mask_array = get_mask_subarray(mask_model, output_model)
    else:
        mask_array = mask_model.dq
    #
    # Bitwise-OR the data pixeldq with the reference file mask
    dq = np.bitwise_or(input_model.pixeldq, mask_array)
    output_model.pixeldq = dq

    output_model.meta.cal_step.dq_init = 'COMPLETE' 

    return output_model

def is_subarray(input_model):

    if (input_model.meta.subarray.xsize==None or
        input_model.meta.subarray.ysize==None):
        raise ValueError('xsize or ysize metadata values not found')

    ncols = input_model.meta.subarray.xsize
    nrows = input_model.meta.subarray.ysize
    instrument = input_model.meta.instrument.name
    if instrument == 'MIRI':
        if ncols < 1032 or nrows < 1024: return True
        return False
    else:
        if ncols < 2048 or nrows < 2048: return True
        return False

def get_mask_subarray(mask_model, output_model):

    if (output_model.meta.subarray.xstart==None or
        output_model.meta.subarray.ystart==None):
        raise ValueError('xstart or ystart metadata values not found')

    xstart = output_model.meta.subarray.xstart
    xstop = xstart + output_model.meta.subarray.xsize - 1
    ystart = output_model.meta.subarray.ystart
    ystop = ystart + output_model.meta.subarray.ysize - 1
    return mask_model.dq[ystart-1:ystop, xstart-1:xstop]

def check_dimensions(input_model):
    #
    # Check that the input model pixeldq attribute has the same dimensions as the
    # image plane of the input model science data
    # If it has dimensions (0,0), create an array of zeros with the same shape as
    # the image plane of the input model
    input_shape = input_model.data.shape
    if input_model.pixeldq.shape != input_shape[-2:]:
        #
        # If the shape is different, then the mask model should have a shape of (0,0)
        # If that's the case, create the array
        if input_model.pixeldq.shape == (0,0):
            input_model.pixeldq = np.zeros((input_shape[-2:])).astype('uint32')
        else:
            print("Pixeldq array has the wrong shape: (%d, %d)" % input_model.pixeldq.shape)
    if input_model.groupdq.shape != input_shape:
        if input_model.groupdq.shape == (0,0,0,0):
            input_model.groupdq = np.zeros((input_shape)).astype('uint8')
        else:
            print("Groupdq array has the wrong shape: (%d, %d, %d, %d)" % input_model.groupdq.shape)
    return
