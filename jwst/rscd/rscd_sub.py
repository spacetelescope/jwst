#
#  Module for applying the RSCD correction to science data
#

import sys
import numpy as np
import logging
from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, rscd_model):
    """
    Short Summary
    -------------
    Applies rscd correction to science arrays, xxx

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd data

    Returns
    -------
    output_model: data model object
        RSCD-corrected science data

    """

    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]       # number of integrations

    sci_ngroups = input_model.data.shape[1]     # number of groups
    # could also grab this information from input_model.meta.exposure.nints (ngroups)



    log.debug("RSCD correction using: nints=%d, ngroups=%d" %
          (sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()

    tau1, scale1, tau2, scale2 = get_rscd_parameters(input_model, rscd_model)
    frame_time = input_model.meta.exposure.frame_time

    # loop over all integrations except the first

    for i in range(1, sci_nints):
        dn_last = get_DNaccumulated_last_int(input_model, i, sci_ngroups)
       # loop over groups in input science data:
        for j in range(sci_ngroups):
           # Apply the correction
            T = (j + 1) * frame_time
            tau = tau1 * frame_time
            eterm = np.exp(-T / tau)
            correction = dn_last * scale1 * eterm
            output.data[i, j] += correction

    return output

def get_rscd_parameters(input_model, rscd_model):
    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name
    print("exposure readpatt,subarray", readpatt, subarray)
    # Read the correct row in table in the reference file to get the coefficients.
    tau1_table = rscd_model.rscd_table.field("tau1")
    scale1_table = rscd_model.rscd_table.field("scale1")
    tau2_table = rscd_model.rscd_table.field("tau2")
    scale2_table = rscd_model.rscd_table.field("scale2")
    readpatt_table = rscd_model.rscd_table.field("readpatt")
    subarray_table = rscd_model.rscd_table.field("subarray")
#    print('tau1',tau1_table)
#    print('readpatt',readpatt_table)
#    print(readpatt_table.size)
#    num_rows = readpatt_table.size


    index = np.asarray(np.where(np.logical_and(readpatt_table == readpatt, subarray_table == subarray)))

    tau1 = tau1_table[index[0]]
    tau2 = tau2_table[index[0]]

    scale1 = scale1_table[index[0]]
    scale2 = scale2_table[index[0]]
    #print('RSCD parameters',tau1,tau2,scale1,scale2)
    return tau1, scale1, tau2, scale2


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
   # print('nrows,ncols',nrows,ncols)
   # print(input_model.data.shape)
    for j in range(nrows):
        for k in range(ncols):
            ramp = input_model.data[i, 0:sci_ngroups - 1, j, k]

            slope, intercept = ols_fit(ramp)
            #print('slope & intercept',slope,intercept)

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
