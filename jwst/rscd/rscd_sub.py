#
#  Module for the RSCD correction for MIRI science data
#

import numpy as np
import logging
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, rscd_model, type):
    """
    Short Summary
    -------------
    if type = baseline the correction sets initial groups in the integration
    to skip
    if type = enhanced the correction corrects the groups using a decay function

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd reference data

    type: string
        type of algorithm ['baseline' or 'enhanced']

    Returns
    -------
    output_model: data model object
        RSCD-corrected science data

    """

    output = input_model.copy()

    # Retrieve the reference parameters for this exposure type
    param = {}
    param = get_rscd_parameters(input_model, rscd_model)

    if param is None:
        log.warning('RSCD correction will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    if type == 'baseline':
        group_skip = param['skip']
        output = correction_skip_groups(input_model, group_skip)
    else:
        log.warning('RSCD correction will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output
        # decay function algorithm updated soon
        #output = correction_decay_function(input_model, param)

    return output


def correction_skip_groups(input_model, nflag):
    """
    Short Summary
    -------------
    Set initial groups in integration to DO_NOT_USE to skip groups
    affected by RSCD effect

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
    sci_ngroups = input_model.data.shape[1]

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # Update the step status, and if ngroups > nflag+3, set all of the GROUPDQ in
    # the first group to 'DO_NOT_USE'
    # the +3 is to ensure there is a slope to be fit including the flagging for
    # the last frame correction
    if sci_ngroups > (nflag + 3):
        output.groupdq[:, 0:nflag :, :] = \
            np.bitwise_or(output.groupdq[:,0:nflag,:,:], dqflags.group['DO_NOT_USE'])
        log.debug(f"RSCD Sub: resetting GROUPDQ in the first {nflag} groups to DO_NOT_USE")
        output.meta.cal_step.rscd = 'COMPLETE'
    else:   # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.rscd = 'SKIPPED'
    return output


def correction_decay_function(input_model, param):
    """
    Short Summary
    -------------
    Applies rscd correction to science arrays using the last group in the last
    integration (after linearity corrected).

    This correction has different correction parameters depending on whether
    the pixel is from an even row or odd row. The first row is defined as an odd
    row. This even/odd row effect is likely a result of the reset electronics
    (MIRI resets in row pairs).

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

    # Create output as a copy of the input science data model
    output = input_model.copy()
    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]       # number of integrations
    sci_ngroups = input_model.data.shape[1]     # number of groups

    # log.debug("RSCD correction using: nints=%d, ngroups=%d" %
    #          (sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()
    return output


def get_rscd_parameters(input_model, rscd_model):


    """
    Read in the parameters from the reference file
    Store the parameters in a param dictionary

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd reference data

    Returns
    -------
    RSCD dictionary of parameters

    """

    # Reference file parameters held in dictionary: param
    param = {}

    # read in the type of data from the input model (FAST,SLOW,FULL,SUBARRAY)
    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name

    # Check for old values of the MIRI LRS slitless subarray name
    # in the science data and change to the new
    if subarray.upper() == 'SUBPRISM': subarray = 'SLITLESSPRISM'

    # read table 1: containing the number of groups to skip
    for tabdata in rscd_model.rscd_group_skip_table:
          subarray_table = tabdata['subarray']
          readpatt_table = tabdata['readpatt']
          group_skip_table = tabdata['group_skip']
          group_skip_table = int(group_skip_table)
          if(subarray_table == subarray and readpatt_table == readpatt):
              param['skip'] = group_skip_table

    # read table 2: General RSCD enhanced parameters
    param['gen'] = {}
    for tabdata in rscd_model.rscd_gen_table:
        readpatt_gen = tabdata['readpatt']
        subarray_gen = tabdata['subarray']
        lower_cutoff_gen = tabdata['lower_cutoff']
        alpha_even_gen = tabdata['alpha_even']
        alpha_odd_gen = tabdata['alpha_even']
        if(subarray_gen == subarray and readpatt_gen == readpatt):
            param['gen']['lower_cutoff'] = lower_cutoff_gen
            param['gen']['lower_alpha_odd'] = alpha_odd_gen
            param['gen']['lower_alpha_even'] = alpha_even_gen

    # read table 3: Enhanced RSCD integration 1 parameters
    param['int1'] = {}
    param['int1']['even'] = {}
    param['int1']['odd'] = {}
    for tabdata in rscd_model.rscd_int1_table:
        readpatt_int1 = tabdata['readpatt']
        subarray_int1= tabdata['subarray']
        rows_int1 = tabdata['rows']
        a0_int1 = tabdata['a0']
        a1_int1 = tabdata['a1']
        a2_int1 = tabdata['a2']
        a3_int1 = tabdata['a3']
        if(subarray_int1 == subarray and readpatt_int1 == readpatt):
            if rows_int1 == 'EVEN':
                param['int1']['even']['a0'] = a0_int1
                param['int1']['even']['a1'] = a1_int1
                param['int1']['even']['a2'] = a2_int1
                param['int1']['even']['a3'] = a3_int1
            if rows_int1 == 'ODD':
                param['int1']['odd']['a0'] = a0_int1
                param['int1']['odd']['a1'] = a1_int1
                param['int1']['odd']['a2'] = a2_int1
                param['int1']['odd']['a3'] = a3_int1

    # read table 4: Enhanced RSCD integration 2 parameters
    param['int2'] = {}
    param['int2']['even'] = {}
    param['int2']['odd'] = {}

    for tabdata in rscd_model.rscd_int2_table:
        readpatt_int2 = tabdata['readpatt']
        subarray_int2= tabdata['subarray']
        rows_int2 = tabdata['rows']
        a0_int2 = tabdata['b0']
        a1_int2 = tabdata['b1']
        a2_int2 = tabdata['b2']
        a3_int2 = tabdata['b3']
        if(subarray_int2 == subarray and readpatt_int2 == readpatt):
             if rows_int2 == 'EVEN':
                 param['int2']['even']['a0'] = a0_int2
                 param['int2']['even']['a1'] = a1_int2
                 param['int2']['even']['a2'] = a2_int2
                 param['int2']['even']['a3'] = a3_int2
             if rows_int2 == 'ODD':
                 param['int2']['odd']['a0'] = a0_int2
                 param['int2']['odd']['a1'] = a1_int2
                 param['int2']['odd']['a2'] = a2_int2
                 param['int2']['odd']['a3'] = a3_int2

    # read table 5: Enhanced RSCD integration 3 parameters
    param['int3'] = {}
    param['int3']['even'] = {}
    param['int3']['odd'] = {}

    for tabdata in rscd_model.rscd_int3_table:
        readpatt_int3 = tabdata['readpatt']
        subarray_int3= tabdata['subarray']
        rows_int3 = tabdata['rows']
        a0_int3 = tabdata['c0']
        a1_int3 = tabdata['c1']
        a2_int3 = tabdata['c2']
        a3_int3 = tabdata['c3']
        if(subarray_int3 == subarray and readpatt_int3 == readpatt):
           if rows_int3 == 'EVEN':
               param['int3']['even']['a0'] = a0_int3
               param['int3']['even']['a1'] = a1_int3
               param['int3']['even']['a2'] = a2_int3
               param['int3']['even']['a3'] = a3_int3
           if rows_int3 == 'ODD':
               param['int3']['odd']['a0'] = a0_int3
               param['int3']['odd']['a1'] = a1_int3
               param['int3']['odd']['a2'] = a2_int3
               param['int3']['odd']['a3'] = a3_int3

    return param
