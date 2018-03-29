#
# Module for handling Data Quality Initialization
#
# This version applies an IRS2-format mask reference file to an IRS2-format
# (3200, 2048) pixeldq array.  Only minor changes were made to the original
# code, in functions is_subarray (use < rather than !=) and check_dimensions
# (pixeldq data type should be uint32).

import logging

import numpy as np

from .. import datamodels
from ..lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# FGS guide star mode exposure types
guider_list = ['FGS_ID-IMAGE', 'FGS_ID-STACK', 'FGS_ACQ1', 'FGS_ACQ2',
               'FGS_TRACK', 'FGS_FINEGUIDE']


def correct_model(input_model, mask_model):
    """DQ Initialize a JWST Model"""
    output_model = do_dqinit(input_model, mask_model)

    return output_model


def do_dqinit(input_model, mask_model):
    """Do the DQ initialization"""

    # Inflate empty DQ array, if necessary
    check_dimensions(input_model)

    # Create output model as copy of input
    output_model = input_model.copy()

    # Extract subarray from reference data, if necessary
    if reffile_utils.ref_matches_sci(output_model, mask_model):
        mask_array = mask_model.dq
    else:
        log.info('Extracting mask subarray to match science data')
        mask_sub_model = reffile_utils.get_subarray_model(output_model,
                                                          mask_model)
        mask_array = mask_sub_model.dq.copy()
        mask_sub_model.close()

    # Set model-specific data quality in output
    if input_model.meta.exposure.type in guider_list:
        dq = np.bitwise_or(input_model.dq, mask_array)
        output_model.dq = dq
    else:
        dq = np.bitwise_or(input_model.pixeldq, mask_array)
        output_model.pixeldq = dq

    output_model.meta.cal_step.dq_init = 'COMPLETE'

    return output_model


def check_dimensions(input_model):
    #
    # Check that the input model pixeldq attribute has the same dimensions as
    # the image plane of the input model science data
    # If it has dimensions (0,0), create an array of zeros with the same shape
    # as the image plane of the input model. For the FGS modes, the
    # GuiderRawModel has only a regular dq array (no pixeldq or groupdq)

    input_shape = input_model.data.shape

    if isinstance(input_model, datamodels.GuiderRawModel):
        if input_model.dq.shape != input_shape[-2:]:

            # If the shape is different, then the mask model should have
            # a shape of (0,0).
            # If that's the case, create the array
            if input_model.dq.shape == (0, 0):
                input_model.dq = np.zeros((input_shape[-2:])).astype('uint32')
            else:
                log.error("DQ array has the wrong shape: (%d, %d)" %
                          input_model.dq.shape)

    else:   # RampModel
        if input_model.pixeldq.shape != input_shape[-2:]:

            # If the shape is different, then the mask model should have
            # a shape of (0,0).
            # If that's the case, create the array
            if input_model.pixeldq.shape == (0, 0):
                input_model.pixeldq = \
                    np.zeros((input_shape[-2:])).astype('uint32')
            else:
                log.error("Pixeldq array has the wrong shape: (%d, %d)" %
                          input_model.pixeldq.shape)

        # Perform the same check for the input model groupdq array
        if input_model.groupdq.shape != input_shape:
            if input_model.groupdq.shape == (0, 0, 0, 0):
                input_model.groupdq = np.zeros((input_shape)).astype('uint8')
            else:
                log.error("Groupdq array has wrong shape: (%d, %d, %d, %d)" %
                          input_model.groupdq.shape)
    return
