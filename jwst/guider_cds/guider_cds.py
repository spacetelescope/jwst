#! /usr/bin/env python
# 
#  guider_cds.py

from __future__ import division
import time
import numpy as np
import logging

from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def guider_cds(model):
    """
    Extended Summary
    ----------------
    Calculate the count rate for each pixel in all data cube sections and all
    integrations.

    For each integration in a given FGS input dataset whose mode is ACQ1, 
    ACQ2, or TRACK, the count rate is the last group minus the first group, 
    divided by the effective integration time.  If the mode is ID, the last 
    group minus the first group is calculated for both integrations; the count 
    rate is then given by the minimum of these two values for each pixel, 
    divided by the effective integration time.  For the FINEGUIDE mode, the 
    count rate is the average of the last 4 groups minus the average of the 
    first 4 groups, divided by the effective integration time.  

    Parameters
    ----------
    model: data model
        input data model, assumed to be of type FGSModel
    """
    # get needed sizes and shapes
    imshape, n_int, effinttm, exp_type = get_dataset_info(model)

    if exp_type[:6] == 'FGS_ID': # force output to have single slice
        new_model = datamodels.GuiderCalModel((1,)+imshape)
    else:
        new_model = datamodels.GuiderCalModel()

    new_model.dq = model.dq
    slope_int_cube = np.zeros(( n_int,) + imshape, dtype = np.float32)

    # loop over data integrations
    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :,:,:]

        if exp_type == 'FGS_FINEGUIDE':
            first_4 = data_sect[:4,:,:].mean(axis=0)
            last_4  = data_sect[-4:,:,:].mean(axis=0)
            slope_int_cube[num_int,:,:] = last_4 - first_4

        elif exp_type[:6] == 'FGS_ID':
            grp_last = data_sect[1, :, :]
            grp_first = data_sect[0, :, :]

            if num_int == 0:
                diff_int0 = grp_last - grp_first
            if num_int == 1:
                diff_int1 = grp_last - grp_first

        else:  # ACQ1, ACQ2, or TRACK
            grp_last = data_sect[1, :, :]
            grp_first = data_sect[0, :, :]
            slope_int_cube[ num_int,:,: ] = grp_last - grp_first
        
    if exp_type[:6] == 'FGS_ID':
        new_model.data[0,:,:] = np.minimum( diff_int1, diff_int0 )/effinttm

    else:  # FINEGUIDE, ACQ1, ACQ2, or TRACK
        new_model.data = slope_int_cube/effinttm

    new_model.update(model)  # ... and add all keys from input

    return new_model


def get_dataset_info(model):
    """
    Short Summary
    -------------
    Extract values for the image shape, the number of integrations, 
    the effective integration time, and the exposure type.

    Parameters
    ----------
    model: instance of Data Model
       DM object for input

    Returns
    -------
    imshape: (int, int) tuple
       shape of 2D image

    n_int: int
       number of integrations

    effinttim: float
       effective integration time

    exp_type: string
        exposure type
    """ 
    instrume = model.meta.instrument.name
    frame_time = model.meta.exposure.frame_time
    ngroups = model.meta.exposure.ngroups

    effinttm = model.meta.exposure.integration_time
    exp_type = model.meta.exposure.type

    n_int = model.data.shape[0]
    nreads = model.data.shape[1]
    asize2 = model.data.shape[2]
    asize1 = model.data.shape[3]

    npix = asize2 * asize1  # number of pixels in 2D array
    imshape = (asize2, asize1)

    log.info('Instrume: %s' % (instrume))
    log.info('Number of integrations: %d' % (n_int))
    log.info('Number of reads: %d' % (nreads))
    log.info('Frame time: %d' % (frame_time))
    log.info('Number of groups per integration: %d' % (ngroups))
    log.info('Effective integration time per group: %s' % (effinttm))
    log.info('Exposure type: %s' % (exp_type))

    return imshape, n_int, effinttm, exp_type
