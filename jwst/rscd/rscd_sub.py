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
    The sole correction is to add the DO_NOT_USE flat to the GROUP data
    quality flags for the first N groups, where N is set by the rscd_model.

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

    if test == 'baseline':
        output = self.set_skip_groups(input)
    else:
        output = self. do_enhanced_correction(input )

    return output


def set_skip_groups(self, input)

    # Save some data params for easy use later
    sci_ngroups = input_model.data.shape[1]

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # get the number of groups to flag
    nflag = get_rscd_parameters(input_model, rscd_model)

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

    # read in the number of groups to skip for simple correction
    readpatt_table = rscd_model.rscd_table['READPATT']
    subarray_table = rscd_model.rscd_table['SUBARRAY']

    # Load the reference table columns of parameters
    readpatt_table = rscd_model.rscd_table['READPATT']
    subarray_table = rscd_model.rscd_table['SUBARRAY']
    rowtype = rscd_model.rscd_table['ROWS']
    tau_table = rscd_model.rscd_table['TAU']
    ascale_table = rscd_model.rscd_table['ASCALE']
    pow_table = rscd_model.rscd_table['POW']
    illum_zp_table = rscd_model.rscd_table['ILLUM_ZP']
    illum_slope_table = rscd_model.rscd_table['ILLUM_SLOPE']
    illum2_table = rscd_model.rscd_table['ILLUM2']
    param3_table = rscd_model.rscd_table['PARAM3']
    crossopt_table = rscd_model.rscd_table['CROSSOPT']
    sat_zp_table = rscd_model.rscd_table['SAT_ZP']
    sat_slope_table = rscd_model.rscd_table['SAT_SLOPE']
    sat_2_table = rscd_model.rscd_table['SAT_2']
    sat_mzp_table = rscd_model.rscd_table['SAT_MZP']
    sat_rowterm_table = rscd_model.rscd_table['SAT_ROWTERM']
    sat_scale_table = rscd_model.rscd_table['SAT_SCALE']

    # Check for old values of the MIRI LRS slitless subarray name
    # in the reference file and change to the new
    for i in range(len(subarray_table)):
        if subarray_table[i].upper() == 'SUBPRISM':
            subarray_table[i] = 'SLITLESSPRISM'

    # Find the matching table row index for even row parameters:
    # the match is based on READPATT, SUBARRAY, and row type (even/odd)
    index = np.asarray(np.where(np.logical_and(readpatt_table == readpatt,
                       np.logical_and(subarray_table == subarray,
                                      rowtype == 'EVEN'))))

    # Check for no match found for EVEN
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=EVEN',
                    readpatt, subarray)
        return None

    # Load the params from the matching table row
    param = {}
    param['even'] = {}
    param['odd'] = {}

    param['even']['tau'] = tau_table[index]
    param['even']['ascale'] = ascale_table[index]
    param['even']['pow'] = pow_table[index]
    param['even']['illum_zp'] = illum_zp_table[index]
    param['even']['illum_slope'] = illum_slope_table[index]
    param['even']['illum2'] = illum2_table[index]
    param['even']['param3'] = param3_table[index]
    param['even']['crossopt'] = crossopt_table[index]
    param['even']['sat_zp'] = sat_zp_table[index]
    param['even']['sat_slope'] = sat_slope_table[index]
    param['even']['sat2'] = sat_2_table[index]
    param['even']['sat_mzp'] = sat_mzp_table[index]
    param['even']['sat_rowterm'] = sat_rowterm_table[index]
    param['even']['sat_scale'] = sat_scale_table[index]

    # Find the matching table row index for ODD row
    index2 = np.asarray(np.where(np.logical_and(readpatt_table == readpatt,
                        np.logical_and(subarray_table == subarray,
                                       rowtype == 'ODD'))))

    # Check for no match found ODD row
    if len(index[0]) == 0:
        log.warning('No matching row found in RSCD reference table for')
        log.warning('READPATT=%s, SUBARRAY=%s, Row type=ODD',
                    readpatt, subarray)
        return None

    # Load the params from the matching table row
    param['odd']['tau'] = tau_table[index2]
    param['odd']['ascale'] = ascale_table[index2]
    param['odd']['pow'] = pow_table[index2]
    param['odd']['illum_zp'] = illum_zp_table[index2]
    param['odd']['illum_slope'] = illum_slope_table[index2]
    param['odd']['illum2'] = illum2_table[index2]
    param['odd']['param3'] = param3_table[index2]
    param['odd']['crossopt'] = crossopt_table[index2]
    param['odd']['sat_zp'] = sat_zp_table[index2]
    param['odd']['sat_slope'] = sat_slope_table[index2]
    param['odd']['sat2'] = sat_2_table[index2]
    param['odd']['sat_mzp'] = sat_mzp_table[index2]
    param['odd']['sat_rowterm'] = sat_rowterm_table[index2]
    param['odd']['sat_scale'] = sat_scale_table[index2]

    return param
    # determine the number of gropus to flag based on the subarray
    #   TBD

    # temp code
    nflag = 3

    return nflag
