#
#  Module for applying the RSCD correction to science data
#

import sys
import numpy as np
import logging
from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def do_correction(input_model, rscd_model):
    """
    Short Summary
    -------------
    Applies rscd correction to science arrays

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd reference data

    Returns
    -------
    output_model: data model object
        RSCD-corrected science data

    """

    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]       # number of integrations
    sci_ngroups = input_model.data.shape[1]     # number of groups
    frame_time = input_model.meta.exposure.frame_time

    log.debug("RSCD correction using: nints=%d, ngroups=%d" %
              (sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # Retrieve the reference parameters for this exposure type
    even, odd = get_rscd_parameters(input_model, rscd_model)

    # Check for valid parameters
    if even is None:
        log.warning('RSCD correction will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    # Unpack the returned parameters
    tau1_even, scale1_even, tau2_even, scale2_even = even
    tau1_odd, scale1_odd, tau2_odd, scale2_odd = odd

    # loop over all integrations except the first
    for i in range(1, sci_nints):
        dn_last = get_DNaccumulated_last_int(input_model, i, sci_ngroups)

        # loop over groups in input science data:
        for j in range(sci_ngroups):

            # Compute the correction factors for even and odd rows
            T = (j + 1) * frame_time
            tau_even = tau1_even * frame_time
            eterm_even = np.exp(-T / tau_even)

            tau_odd = tau1_odd * frame_time
            eterm_odd = np.exp(-T / tau_odd)

            # Apply the corrections to even and odd rows:
            # the first row is defined as odd (python index 0)
            # the second row is the first even row (python index of 1)
            correction_odd = dn_last[0::2, :] * scale1_odd * eterm_odd
            correction_even = dn_last[1::2, :] * scale1_even * eterm_even

            output.data[i, j, 0::2, :] += correction_odd
            output.data[i, j, 1::2, :] += correction_even

    output.meta.cal_step.rscd = 'COMPLETE'

    return output


def get_rscd_parameters(input_model, rscd_model):

    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name

    # Load the reference table columns of parameters
    tau1_table = rscd_model.rscd_table['TAU1']
    scale1_table = rscd_model.rscd_table['SCALE1']
    tau2_table = rscd_model.rscd_table['TAU2']
    scale2_table = rscd_model.rscd_table['SCALE2']
    readpatt_table = rscd_model.rscd_table['READPATT']
    subarray_table = rscd_model.rscd_table['SUBARRAY']
    rowtype = rscd_model.rscd_table['rows']

    # Find the matching table row index for even row parameters:
    # the match is based on READPATT, SUBARRAY, and row type (even/odd)
    index = np.asarray(np.where(np.logical_and(readpatt_table == readpatt,
                       np.logical_and(subarray_table == subarray, rowtype == 'EVEN'))))

    # Check for no match found
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=EVEN',
                    readpatt, subarray)
        return None, None

    # Load the params from the matching table row
    tau1_even = tau1_table[index]
    tau2_even = tau2_table[index]
    scale1_even = scale1_table[index]
    scale2_even = scale2_table[index]

    # Find the matching table row index for odd row parameters:
    index2 = np.asarray(np.where(np.logical_and(readpatt_table == readpatt,
                        np.logical_and(subarray_table == subarray, rowtype == 'ODD'))))

    # Check for no match found
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=ODD',
                    readpatt, subarray)
        return None, None

    # Load the params from the matching table row
    tau1_odd = tau1_table[index2]
    tau2_odd = tau2_table[index2]
    scale1_odd = scale1_table[index2]
    scale2_odd = scale2_table[index2]

    even = tau1_even, scale1_even, tau2_even, scale2_even
    odd = tau1_odd, scale1_odd, tau2_odd, scale2_odd

    return even, odd


def get_DNaccumulated_last_int(input_model, i, sci_ngroups):

    # Find the accumulated DN from the last integration
    # need to add skipping N frames (frames affected by reset)
    # Add check to make sure we have enough frames left to do a fit

    nrows = input_model.data.shape[2]
    ncols = input_model.data.shape[3]

    # last frame affected by "last frame" effect - use second to last frame
    # may want to extrapolate to last frame
    # we may want to check if data has saturated
    dn_lastframe = input_model.data[i - 1][sci_ngroups - 2]

    dn_accumulated = dn_lastframe.copy() * 0.0
    for j in range(nrows):
        for k in range(ncols):
            ramp = input_model.data[i, 0:sci_ngroups - 1, j, k]
            slope, intercept = ols_fit(ramp)
            dn_accumulated[j, k] = dn_lastframe[j, k] - intercept

    return dn_accumulated


def ols_fit(y):

    shape = y.shape
    nelem = float(len(y))

    x = np.arange(shape[0], dtype=np.float64)
    xshape = list(shape)
    for i in range(1, len(shape)):
        xshape[i] = 1
    x = x.reshape(xshape)

    mean_y = y.mean(axis=0)
    mean_x = x[-1] / 2.
    sum_x2 = (x**2).sum(axis=0)
    sum_xy = (x * y).sum(axis=0)
    slope = (sum_xy - nelem * mean_x * mean_y) / \
            (sum_x2 - nelem * mean_x**2)
    intercept = mean_y - slope * mean_x

    return (slope, intercept)
