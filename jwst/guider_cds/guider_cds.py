#! /usr/bin/env python
#
#  guider_cds.py

import numpy as np
import logging

from .. import datamodels
from ..lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# When an appropriate ref file for the gain (read noise) can not be retrieved
#  for the calculation of the variances, and then the ERR array constant
#  value equal to the mean of the SCI values in this accessible ref file will
#  be used:
DEFAULT_GAIN = 1.9713863
DEFAULT_READNOISE = 21.0
DEFAULT_GAIN_FILE = "jwst_fgs_gain_0010.fits"
DEFAULT_READNOISE_FILE = "jwst_fgs_readnoise_0000.fits"


def guider_cds(model, gain_model, readnoise_model):
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
    divided by the group time.  The variances are calculated using gain and
    read noise values obtained from reference files if possible, otherwise,
    default representative constant values are used. From the variances the ERR
    array is calculated.

    Parameters
    ----------
    model: data model
        input data model, assumed to be of type GuiderRawModel

    gain_model: instance of gain Model or None
        gain for all pixels

    readnoise_model: instance of readnoise Model or None
        readnoise for all pixels
    """
    # get needed sizes and shapes
    imshape, n_int, grp_time, exp_type = get_dataset_info(model)

    # get gain and readnoise arrays to calculate ERR array
    gain_arr, readnoise_arr = get_ref_arr(model, imshape, gain_model,
                                          readnoise_model)

    if exp_type[:6] == 'FGS_ID':  # force output to have single slice
        new_model = datamodels.GuiderCalModel((1,) + imshape)
    else:
        new_model = datamodels.GuiderCalModel()

    # set up output data arrays
    slope_int_cube = np.zeros((n_int,) + imshape, dtype=np.float32)
    if exp_type[:6] == 'FGS_ID':
        var_rn = np.zeros((1,) + imshape, dtype=np.float32)
        var_pn = np.zeros((1,) + imshape, dtype=np.float32)
    else:
        var_rn = np.zeros((n_int,) + imshape, dtype=np.float32)
        var_pn = np.zeros((n_int,) + imshape, dtype=np.float32)
    new_model.dq = model.dq

    # loop over data integrations
    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]
        if exp_type == 'FGS_FINEGUIDE':
            first_4 = data_sect[:4, :, :].mean(axis=0)
            last_4 = data_sect[-4:, :, :].mean(axis=0)
            slope_int_cube[num_int, :, :] = last_4 - first_4

            var_rn[num_int, :, :] = 2 * (readnoise_arr / grp_time)**2
            var_pn[num_int, :, :] = slope_int_cube[num_int, :, :] / (gain_arr * grp_time)

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
            var_rn[num_int, :, :] = 2 * (readnoise_arr / grp_time)**2
            var_pn[num_int, :, :] = slope_int_cube[num_int, :, :] / (gain_arr * grp_time)

    if exp_type[:6] == 'FGS_ID':
        new_model.data[0, :, :] = np.minimum(diff_int1, diff_int0) / grp_time
    else:  # FINEGUIDE, ACQ1, ACQ2, or TRACK
        new_model.data = slope_int_cube / grp_time

    # set err to sum of variances in quadrature
    new_model.err = (var_rn * var_rn + var_pn * var_pn) ** 0.5

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


def get_ref_arr(model, imshape, gain_model, readnoise_model):
    """
    Short Summary
    -------------
    Create gain and readnoise files from their reference files

    Parameters
    ----------
    model: data model
        input data model, assumed to be of type GuiderRawModel

    imshape: (int, int) tuple
       shape of 2D image

    gain_model: instance of gain Model or None
        gain for all pixels

    readnoise_model: instance of readnoise Model or None
        readnoise for all pixels

    Returns
    -------
    gain_arr: ndarray, 2-D, float
        gain values from the ref file, or default value if ref file unavailable

    readnoise_arr: ndarray, 2-D, float
        readnoise values from the ref file, or default value if ref file
        unavailable
    """
    # if no gain_model was retrieved, use the default value for all pixels
    if gain_model is None:  # no ref file, so set all pixel's gain to constant default
        gain_arr = np.zeros(imshape, dtype=np.float32) + DEFAULT_GAIN

        log.info('The appropriate GAIN reference file could not be retrieved,')
        log.info('so instead, from the reference file %s', DEFAULT_GAIN_FILE)
        log.info('all gain values will be set to the mean of the SCI values (%s)', DEFAULT_GAIN)
    else:
        # extract subarray from gain reference file, if necessary
        if reffile_utils.ref_matches_sci(model, gain_model):
            gain_arr = gain_model.data
        else:
            log.info('Extracting reference file subarray to match science data')
            ref_sub_model = reffile_utils.get_subarray_model(model, gain_model)
            gain_arr = ref_sub_model.data.copy()
            ref_sub_model.close()

    # if no readnoise_model was retrieved, use the default value for all pixels
    if readnoise_model is None:  # no ref file, so set all pixel's readnoise to constant
        readnoise_arr = np.zeros(imshape, dtype=np.float32) + DEFAULT_READNOISE

        log.info('The appropriate READNOISE reference file could not be retrieved,')
        log.info('so instead, from the reference file %s', DEFAULT_READNOISE_FILE)
        log.info('all read noise values will be set to the mean of the SCI values (%s)', DEFAULT_READNOISE)
    else:
        # extract subarray from readnoise reference file, if necessary
        if reffile_utils.ref_matches_sci(model, readnoise_model):
            readnoise_arr = readnoise_model.data
        else:
            log.info('Extracting readnoise reference file subarray to match science data')
            ref_sub_model = reffile_utils.get_subarray_model(model, readnoise_model)
            readnoise_arr = ref_sub_model.data.copy()
            ref_sub_model.close()

    return gain_arr, readnoise_arr


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
