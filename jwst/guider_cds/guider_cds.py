#! /usr/bin/env python
#
#  guider_cds.py

import numpy as np
import logging

from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def guider_cds(model):
    """
    Extended Summary
    ----------------
    Calculate the count rate for each pixel in all integrations.

    For each integration in a given FGS guider dataset whose mode is ACQ1,
    ACQ2, or TRACK, the count rate is the last group minus the first group,
    divided by the effective integration time.  If the mode is ID, the last
    group minus the first group is calculated for both integrations; the count
    rate is then given by the minimum of these two values for each pixel,
    divided by the group time.  For the FINEGUIDE mode, the count rate is the
    average of the last 4 groups minus the average of the first 4 groups,
    divided by the group time.

    Parameters
    ----------
    model: data model
        input data model, assumed to be of type GuiderRawModel
    """

    # get needed sizes and shapes
    imshape, n_int, grp_time, exp_type = get_dataset_info(model)

    if exp_type[:6] == 'FGS_ID':  # force output to have single slice
        new_model = datamodels.GuiderCalModel((1,) + imshape)
    else:
        new_model = datamodels.GuiderCalModel()

    # set up output data arrays
    slope_int_cube = np.zeros((n_int,) + imshape, dtype=np.float32)
    if exp_type[:6] == 'FGS_ID':
        err = np.zeros((1,) + imshape, dtype=np.float32)
    else:
        err = np.zeros((n_int,) + imshape, dtype=np.float32)
    new_model.err = err
    new_model.dq = model.dq

    # loop over data integrations
    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]

        if exp_type == 'FGS_FINEGUIDE':
            first_4 = data_sect[:4, :, :].mean(axis=0)
            last_4 = data_sect[-4:, :, :].mean(axis=0)
            slope_int_cube[num_int, :, :] = last_4 - first_4

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
            slope_int_cube[num_int, :, :] = grp_last - grp_first

    if exp_type[:6] == 'FGS_ID':
        new_model.data[0, :, :] = np.minimum(diff_int1, diff_int0) / grp_time
    else:  # FINEGUIDE, ACQ1, ACQ2, or TRACK
        new_model.data = slope_int_cube / grp_time

    # Add all table extensions to be carried over to output
    if len(model.planned_star_table):
        new_model.planned_star_table = model.planned_star_table
    if len(model.flight_star_table):
        new_model.flight_star_table = model.flight_star_table
    if len(model.pointing_table):
        new_model.pointing_table = model.pointing_table
    if len(model.centroid_table):
        new_model.centroid_table = model.centroid_table
    if len(model.track_sub_table):
        new_model.track_sub_table = model.track_sub_table

    # copy all meta data from input to output model
    new_model.update(model)

    # Update BUNIT to reflect count rate
    new_model.meta.bunit_data = 'DN/s'

    return new_model


def get_dataset_info(model):
    """
    Short Summary
    -------------
    Extract values for the image shape, the number of integrations,
    the group time, and the exposure type.

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

    grp_time: float
       group time

    exp_type: string
        exposure type
    """
    instrume = model.meta.instrument.name
    frame_time = model.meta.exposure.frame_time
    ngroups = model.meta.exposure.ngroups
    grp_time = model.meta.exposure.group_time
    exp_type = model.meta.exposure.type

    n_int = model.data.shape[0]
    asize2 = model.data.shape[2]
    asize1 = model.data.shape[3]

    imshape = (asize2, asize1)

    log.info('Instrument: %s' % (instrume))
    log.info('Exposure type: %s' % (exp_type))
    log.info('Number of integrations: %d' % (n_int))
    log.info('Number of groups per integration: %d' % (ngroups))
    log.info('Group time: %s' % (grp_time))
    log.info('Frame time: %10.5f' % (frame_time))

    return imshape, n_int, grp_time, exp_type
