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
    The last frame value from the previous integration is calculated two ways:
    1. using second and third to last frames to extrapolated to the last frame
    2. using the non saturating data, fit the data and expolate to last frame
    Because of the uncertainity of how well effects in th early part of the integration
    are corrected in the previous integration (reset anomaly, rscd effects, persistence)
    the lastframe determined from the second and third to last frames is considered
    a better estimate than that derived from a fit to the ramp.
    The last frame derived from fitting the non-saturating data is used in the correction
    if the previous integration saturated. This fit is expolated past saturation to estimate
    what the total number of electrons would of been collected was. 

    This correction has different correction parameters depending on if the pixel is from an
    even row or odd row. The first row is define as an odd row. This even/odd row effect is
    likely a result of the reset electronics (MIRI resets in row pairs).

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

    # Check for valid parameters
    if sci_ngroups < 2:
        log.warning('RSCD correction will be skipped, only 1 Group need at least 2')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    # Retrieve the reference parameters for this exposure type
    param = get_rscd_parameters(input_model, rscd_model)

    if param is None:
        log.warning('RSCD correction will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    # Determine the parameters that only rely on ngroups
    
    ngroups2 = sci_ngroups * sci_ngroups
    b1_even = param['even']['ascale'] * (param['even']['illum_zp'] +
                                         param['even']['illum_slope']*sci_ngroups +
                                         param['even']['illum2']*ngroups2)

    b1_odd = param['odd']['ascale'] * (param['odd']['illum_zp'] +
                                         param['odd']['illum_slope']*sci_ngroups +
                                         param['odd']['illum2']*ngroups2)

    sat_final_slope_even = (param['even']['sat_zp'] + param['even']['sat_slope']*sci_ngroups +
                            param['even']['sat2']*ngroups2 + param['even']['sat_rowterm'])

    sat_final_slope_odd = (param['odd']['sat_zp'] + param['odd']['sat_slope']*sci_ngroups +
                            param['odd']['sat2']*ngroups2 + param['odd']['sat_rowterm'])


    b2_even = np.asscalar(param['even']['pow'])   
    b2_odd = np.asscalar(param['odd']['pow'])
    b3_even = np.asscalar(param['even']['param3'])
    b3_odd = np.asscalar(param['odd']['param3'])
    crossopt_even = np.asscalar(param['even']['crossopt'])
    crossopt_odd = np.asscalar(param['odd']['crossopt'])
    sat_mzp_even = np.asscalar(param['even']['sat_mzp'])
    sat_mzp_odd = np.asscalar(param['odd']['sat_mzp'])
    sat_scale_even = np.asscalar(param['even']['sat_scale'])
    sat_scale_odd = np.asscalar(param['odd']['sat_scale'])
    tau_even = np.asscalar(param['even']['tau'])
    tau_odd = np.asscalar(param['odd']['tau'])
    # loop over all integrations except the first
    for i in range(1, sci_nints):

        sat,is_ref,dn_last23,dn_lastfit = get_DNaccumulated_last_int(input_model, i, sci_ngroups)
        lastframe_even = dn_last23[1::2,:]
        lastframe_odd = dn_last23[0::2,:]

        correction_even = lastframe_even.copy() * 0.0
        correction_odd = lastframe_odd.copy ()* 0.0
        factor2_even = lastframe_even.copy ()* 0.0
        factor2_odd = lastframe_odd.copy ()* 0.0
        a1_even = lastframe_even.copy ()* 0.0
        a1_odd = lastframe_odd.copy ()* 0.0

        counts2_even = lastframe_even - crossopt_even
        counts2_odd =  lastframe_odd - crossopt_odd
        
        counts2_even[np.where(counts2_even< 0)] = 0.0
        counts2_odd[ np.where(counts2_odd< 0)] = 0.0

        # Find where counts2 > 0 and is finite
        good_even = np.where(counts2_even >  0 & np.isfinite(counts2_even))
        good_odd = np.where(counts2_odd  > 0 & np.isfinite(counts2_odd))
        #______________________________________________________________________
        # even row values
        factor2_even[good_even] = 1.0/(np.exp(counts2_even[good_even]/b3_even) -1)
        a1_even = b1_even *(np.power(counts2_even,b2_even)) * factor2_even
        #______________________________________________________________________
        # odd row values
        factor2_odd[good_odd] = 1.0/(np.exp(counts2_odd[good_odd]/b3_odd) -1)
        a1_odd = b1_odd *(np.power(counts2_odd,b2_odd)) * factor2_odd

        #______________________________________________________________________
        # SATURATED DATA 
        counts3_even = dn_lastfit[1::2,:] * sat_scale_even                    
        counts3_odd = dn_lastfit[0::2,:] * sat_scale_odd

        a1_sat_even =sat_final_slope_even*counts3_even + sat_mzp_even
        a1_sat_odd =sat_final_slope_odd*counts3_odd + sat_mzp_odd

        sat_even = sat[1::2,:]
        sat_odd = sat[0::2,:]

        # loop over groups in input science data:
        for j in range(sci_ngroups):
                                
            # Compute the correction factors for even and odd rows
            T = (j + 1) 
                            
            eterm_even = np.exp(-T / tau_even)
            eterm_odd = np.exp(-T / tau_odd)

            # Apply the corrections to even and odd rows:
            # the first row is defined as odd (python index 0)
            # the second row is the first even row (python index of 1)
            correction_odd = lastframe_odd * a1_odd * 0.01*eterm_odd
            correction_even = lastframe_even * a1_even * 0.01* eterm_even

            correction_sat_odd = lastframe_odd *  a1_sat_odd * 0.01*  eterm_odd
            correction_sat_even = lastframe_even * a1_sat_even* 0.01 * eterm_even
            sat_index_even  = np.where(sat_even)
            sat_index_odd  = np.where(sat_odd)
            correction_even[sat_index_even] = correction_sat_even[sat_index_even]
            correction_odd[sat_index_odd] = correction_sat_odd[sat_index_odd]

            output.data[i, j, 0::2, :] += correction_odd
            output.data[i, j, 1::2, :] += correction_even


    output.meta.cal_step.rscd = 'COMPLETE'

    return output


def get_rscd_parameters(input_model, rscd_model):
    
    """
    Read in the parameters from the reference file
    Seperate these parameters based on even and odd rows
    
    Store the parameters in a param dictionary
    """
    # read in the type of data from the input model (FAST,SLOW,FULL,SUBARRAY)
    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name

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
    param = {}
    param['even'] = {}
    param['odd'] = {}
    
    param['even']['tau'] = tau_table[index]
    param['even']['ascale'] = ascale_table[index]
    param['even']['pow']  = pow_table[index]
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
    param['odd']['tau'] = tau_table[index2]
    param['odd']['ascale'] = ascale_table[index2]
    param['odd']['pow']  = pow_table[index2]
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


def get_DNaccumulated_last_int(input_model, i, sci_ngroups):

    """
    Find the accumulated DN from the last integration
    This data should already have the Reset Anomaly correction
    applied (if not - should we skip frames at the beginning ?) 

    a Check has already been made to make sure we have at least 
    4 frames 

    Parameters
    ----------
    input_model: input ramp data
    i: integration #
    sci_ngroups: number of frames/integration
    
    return values
    -------------
    sat: is the previous integration for this pixel saturated: yes/no
    is_ref: boolean if pixel is a reference pixel. Do not apply correction to these pixels
            we set lastframe = 0 for these pixels 
    dn_lastframe_23: the extrapolated last frame using second and third to last frames
    dn_lastfrane_fit: the extrapolated last frame using the fit to the entire ramp
   
    """
 
    nrows = input_model.data.shape[2]
    ncols = input_model.data.shape[3]

    dn_lastframe2 = input_model.data[i - 1][sci_ngroups - 2]
    dn_lastframe3 = input_model.data[i - 1][sci_ngroups - 3]
    dn_lastframe23 = dn_lastframe2.copy() * 0.0
    dn_lastframe_fit = dn_lastframe2.copy() * 0.0
    saturated = np.full((nrows,ncols),False)
    is_ref = saturated.copy()

    diff = dn_lastframe2 - dn_lastframe3
    dn_lastframe23 = dn_lastframe2 + diff
    # check if pixel is saturated
    saturated_flag = datamodels.dqflags.group['SATURATED']
    ref_flag = datamodels.dqflags.pixel['REFERENCE_PIXEL']
#TODO make this section more efficient by not looping over each
# pixel
    for j in range(nrows):
        for k in range(ncols):
            pixeldq = input_model.pixeldq[j,k]
            ref = np.bitwise_and(pixeldq,ref_flag) == 0
            if ref == True:
                is_ref[j,k] = True
                dn_lastframe23[j,k] = 0.0

            else : # pixel is not a reference pixels. Make a correction
                is_ref[j,k] = False
                # starting on second frame and going to second to last frame
                # pull out the pixel ramp. The first and last frame are
                # heavily effected by detector effects and may not be
                # adequately corrected for this simple ols fit.

                ramp = input_model.data[i-1, 1:sci_ngroups - 1, j, k]
                groupdq = input_model.groupdq[i-1,1:sci_ngroups-1,j,k]
                satmask = (groupdq == saturated_flag)
                yessat = satmask.any()
                saturated[j,k] = yessat
            
                slope, intercept,ngood = ols_fit(ramp,groupdq)
                if slope !=0: 
                    dn_lastframe_fit[j,k] = slope*sci_ngroups + intercept
                else: 
                    dn_lastframe_fit[j,k] = dn_lastframe23[j,k]

    return saturated, is_ref,dn_lastframe23,dn_lastframe_fit


def ols_fit(y,dq):

    """
    An estimation of the lastframe Value from the previous integration is needed
    for the RSCD correction.
    This routine does a simple ordinary least squares fit to non-saturating data
    """

    saturated_flag = datamodels.dqflags.group['SATURATED']
    shape = y.shape
    ngood = 0

    x = np.arange(shape[0], dtype=np.float64)
    xshape = list(shape)

    for i in range(1, len(shape)):
        xshape[i] = 1
    x = x.reshape(xshape)

    good_data = np.where(np.bitwise_and(dq,saturated_flag)==0)
    ngood = len(good_data[0])
    slope = 0.0
    intercept = 0.0 

    if ngood >= 3:
        xuse = x[good_data]
        yuse = y[good_data]
        nelem = float(len(yuse))
        mean_y = yuse.mean(axis=0)
        mean_x = xuse[-1] / 2.
        sum_x2 = (xuse**2).sum(axis=0)
        sum_xy = (xuse * yuse).sum(axis=0)
        slope = (sum_xy - nelem * mean_x * mean_y) / \
            (sum_x2 - nelem * mean_x**2)
        intercept = mean_y - slope * mean_x
        
    return (slope, intercept,ngood)
