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

    # Retrieve the reference parameters for this exposure type
    param = get_rscd_parameters(input_model, rscd_model)

    if param is None:
        log.warning('RSCD correction will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    output = input_model.copy()
    if type == 'baseline':
        group_skip = param['skip'][0]
        output = self.correction_skip_groups(input_model,group_skip)
    else:
        output = self.correction_decay_function(input_model,param)

    return output


def correction_skip_groups(input_model,nflag):

    # Save some data params for easy use later
    sci_ngroups = input_model.data.shape[1]

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # Update the step status, and if ngroups > nflag+3, set all of the GROUPDQ in
    # the first group to 'DO_NOT_USE'
    # the +3 is to ensure there is a slope to be fit including the flagging for
    # the last frame correction
    if sci_ngroups > (nflag + 3):
        output.groupdq[:, 0, :, :] = \
            np.bitwise_or(output.groupdq[:,0:nflag,:,:], dqflags.group['DO_NOT_USE'])
        log.debug(f"RSCD Sub: resetting GROUPDQ in the first {nflag} groups to DO_NOT_USE")
        output.meta.cal_step.rscd = 'COMPLETE'
    else:   # too few groups
        log.warning("Too few groups to apply correction")
        log.warning("Step will be skipped")
        output.meta.cal_step.rscd = 'SKIPPED'


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
    RSCD tables of parameters

    """

    # read in the type of data from the input model (FAST,SLOW,FULL,SUBARRAY)
    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name

    # Check for old values of the MIRI LRS slitless subarray name
    # in the science data and change to the new
    if subarray.upper() == 'SUBPRISM': subarray = 'SLITLESSPRISM'

    # read table 1: containing the number of groups to skip
    subarray_table1 = rscd_model.rscd_group_skip_table['subarray']
    readpatt_table1 = rscd_model.rscd_group_skip_table['readpatt']
    group_skip_table1 = rscd_model.rscd_group_skip_table['group_skip']


    # read table 2: General RSCD enahnced parameters
    readpatt_gen = rscd_model.rscd_gen_table['readpatt']
    subarray_gen = rscd_model.rscd_gen_table['subarray']
    lower_cutoff_gen = rscd_model.rscd_gen_table['lower_cutoff']
    alpha_even_gen = rscd_model.rscd_gen_table['alpha_even']
    alpha_odd_gen = rscd_model.rscd_gen_table['alpha_even']

    # read table 3: General RSCD integration 1 parameters
    readpatt_int1 = rscd_model.rscd_int1_table['readpatt']
    subarray_int1= rscd_model.rscd_int1_table['subarray']
    rows_int1 = rscd_model.rscd_int1_table['rows']
    a0_int1 = rscd_model.rscd_int1_table['a0']
    a1_int1 = rscd_model.rscd_int1_table['a1']
    a2_int1 = rscd_model.rscd_int1_table['a2']
    a3_int1 = rscd_model.rscd_int1_table['a3']

    # read table 3: General RSCD integration 2 parameters
    readpatt_int2 = rscd_model.rscd_int2_table['readpatt']
    subarray_int2= rscd_model.rscd_int2_table['subarray']
    rows_int2 = rscd_model.rscd_int2_table['rows']
    a0_int2 = rscd_model.rscd_int2_table['b0']
    a1_int2 = rscd_model.rscd_int2_table['b1']
    a2_int2 = rscd_model.rscd_int2_table['b2']
    a3_int2 = rscd_model.rscd_int2_table['b3']

    # read table 3: General RSCD integration 1 parameters
    readpatt_int3 = rscd_model.rscd_int3ç_table['readpatt']
    subarray_int3= rscd_model.rscd_int3_table['subarray']
    rows_int3 = rscd_model.rscd_int3_table['rows']
    a0_int3 = rscd_model.rscd_int3_table['c0']
    a1_int3 = rscd_model.rscd_int3_table['c1']
    a2_int3 = rscd_model.rscd_int3_table['c2']
    a3_int3 = rscd_model.rscd_int3_table['c3']

    # Check for old values of the MIRI LRS slitless subarray name
    # in the reference file and change to the new
 #   for i in range(len(subarray_table)):
 #       if subarray_table[i].upper() == 'SUBPRISM':
 #           subarray_table[i] = 'SLITLESSPRISM'

    # for the Group Skip match data to table with Readpatt and subarray
    index = np.asarray(np.where(np.logical_and(readpatt_table1 == readpatt,
                       subarray_table1 == subarray)))

    # Check for no match found for readpatt and subarray
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s', readpatt, subarray)
        return None

     #Load the params from the matching table row
    param = {}
    param['skip'] = {}
    param['skip'] = group_skip_table1[index]

     # fill in the RSCD enhanced gen values 
    # for the RSCD Gen  match data to table with Readpatt and subarray
    index = np.asarray(np.where(np.logical_and(readpatt_gen == readpatt,
                       subarray_gen == subarray)))

    param['gen'] = {}
    param['gen']['lower_cutoff'] = lower_cutoff_gen[index]
    param['gen']['lower_alpha_odd'] = alpha_odd_gen[index]
    param['gen']['lower_alpha_even'] = alpha_even_gen[index]

     # fill in the RSCD enhanced integration 1  
    # for the RSCD Gen  match data to table with Readpatt, subarray, row

    index = np.where(np.logical_and(readpatt_int1 == readpatt,
                                     np.logical_and(subarray_int == subarray,
                                                    rows_int == 'EVEN')))
    # Check for no match found for EVEN
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=EVEN',
                    readpatt, subarray)
        return None

    param['int1'] = {}
    param['int1']['even'] = {}
    param['int1']['even']['a0'] = a0_int1[index]
    param['int1']['even']['a1'] = a1_int1[index]
    param['int1']['even']['a2'] = a2_int1[index]
    param['int1']['even']['a3'] = a3_int1[index]

    index = np.where(np.logical_and(readpatt_int1 == readpatt,
                                     np.logical_and(subarray_int1 == subarray,
                                                    rows_int1 == 'ODD')))

    # Check for no match found for ODD
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=ODD',
                    readpatt, subarray)
        return None
    param['int1'] = {}
    param['int1']['odd'] = {}
    param['int1']['odd']['a0'] = a0_int1[index]
    param['int1']['odd']['a1'] = a1_int1[index]
    param['int1']['odd']['a2'] = a2_int1[index]
    param['int1']['odd']['a3'] = a3_int1[index]

     # fill in the RSCD enhanced integration 2  
    # for the RSCD Gen  match data to table with Readpatt, subarray, row

    index = np.where(np.logical_and(readpatt_int2 == readpatt,
                                     np.logical_and(subarray_int2 == subarray,
                                                    rows_int2 == 'EVEN')))

    # Check for no match found for EVEN
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=EVEN',
                    readpatt, subarray)
        return None
    param['int2'] = {}
    param['int2']['even'] = {}
    param['int2']['even']['a0'] = a0_int2[index]
    param['int2']['even']['a1'] = a1_int2[index]
    param['int2']['even']['a2'] = a2_int2[index]
    param['int2']['even']['a3'] = a3_int2[index]

    index = np.where(np.logical_and(readpatt_int2 == readpatt,
                                     np.logical_and(subarray_int2 == subarray,
                                                    rows_int2 == 'ODD')))

    # Check for no match found for ODD
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=ODD',
                    readpatt, subarray)
        return None
    param['int2'] = {}
    param['int2']['odd'] = {}
    param['int2']['odd']['a0'] = a0_int2[index]
    param['int2']['odd']['a1'] = a1_int2[index]
    param['int2']['odd']['a2'] = a2_int2[index]
    param['int2']['odd']['a3'] = a3_int2[index]

     # fill in the RSCD enhanced integration 3
    # for the RSCD Gen  match data to table with Readpatt, subarray, row

    index = np.where(np.logical_and(readpatt_int3 == readpatt,
                                     np.logical_and(subarray_int3 == subarray,
                                                    rows_int3 == 'EVEN')))
    # Check for no match found for EVEN
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=EVEN',
                    readpatt, subarray)
        return None
    param['int3'] = {}
    param['int3']['even'] = {}
    param['int3']['even']['a0'] = a0_int3[index]
    param['int3']['even']['a1'] = a1_int3[index]
    param['int3']['even']['a2'] = a2_int3[index]
    param['int3']['even']['a3'] = a3_int3[index]

    index = np.where(np.logical_and(readpatt_int3 == readpatt,
                                     np.logical_and(subarray_int3 == subarray,
                                                    rows_int3 == 'ODD')))
    # Check for no match found for ODD
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=ODD',
                    readpatt, subarray)
        return None
    param['int3'] = {}
    param['int3']['odd'] = {}
    param['int3']['odd']['a0'] = a0_int3[index]
    param['int3']['odd']['a1'] = a1_int3[index]
    param['int3']['odd']['a2'] = a2_int3[index]
    param['int3']['odd']['a3'] = a3_int3[index]


    return param
