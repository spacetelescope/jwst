#
#  Module for applying the RSCD correction to science data
#

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
    Because of the uncertainity of how well effects in th early part of the
    integration are corrected in the previous integration (reset anomaly, rscd
    effects, persistence) the lastframe determined from the second and third to
    last frames is considered a better estimate than that derived from a fit to
    the ramp.
    The last frame derived from fitting the non-saturating data is used in the
    correction if the previous integration saturated. This fit is expolated
    past saturation to estimate what the total number of electrons would have
    been.

    This correction has different correction parameters depending on whether
    the pixel is from an even row or odd row. The first row is define as an odd
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

    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]       # number of integrations
    sci_ngroups = input_model.data.shape[1]     # number of groups

    log.debug("RSCD correction using: nints=%d, ngroups=%d" %
              (sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # Check for valid parameters
    if sci_ngroups < 2:
        log.warning('RSCD correction requires > 1 group per integration')
        log.warning('Step will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    # Retrieve the reference parameters for this exposure type
    param = get_rscd_parameters(input_model, rscd_model)

    if param is None:
        log.warning('RSCD correction will be skipped')
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    # Determine the parameters that rely only on ngroups
    ngroups2 = sci_ngroups * sci_ngroups
    b1_even = param['even']['ascale'] * (
        param['even']['illum_zp'] +
        param['even']['illum_slope'] * sci_ngroups +
        param['even']['illum2'] * ngroups2)

    b1_odd = param['odd']['ascale'] * (
        param['odd']['illum_zp'] +
        param['odd']['illum_slope'] * sci_ngroups +
        param['odd']['illum2'] * ngroups2)

    sat_final_slope_even = (
        param['even']['sat_zp'] + param['even']['sat_slope'] * sci_ngroups +
        param['even']['sat2'] * ngroups2 + param['even']['sat_rowterm'])

    sat_final_slope_odd = (
        param['odd']['sat_zp'] + param['odd']['sat_slope'] * sci_ngroups +
        param['odd']['sat2'] * ngroups2 + param['odd']['sat_rowterm'])

    b2_even = param['even']['pow'].item()
    b2_odd = param['odd']['pow'].item()
    b3_even = param['even']['param3'].item()
    b3_odd = param['odd']['param3'].item()
    crossopt_even = param['even']['crossopt'].item()
    crossopt_odd = param['odd']['crossopt'].item()
    sat_mzp_even = param['even']['sat_mzp'].item()
    sat_mzp_odd = param['odd']['sat_mzp'].item()
    sat_scale_even = param['even']['sat_scale'].item()
    sat_scale_odd = param['odd']['sat_scale'].item()
    tau_even = param['even']['tau'].item()
    tau_odd = param['odd']['tau'].item()

    # loop over all integrations except the first
    mdelta = int(sci_nints / 10) + 1
    for i in range(1, sci_nints):
        if ((i + 1) % mdelta) == 0:
            log.info(' Working on integration %d', i + 1)

        sat, dn_last23, dn_lastfit = \
            get_DNaccumulated_last_int(input_model, i, sci_ngroups)

        lastframe_even = dn_last23[1::2, :]
        lastframe_odd = dn_last23[0::2, :]

        correction_even = lastframe_even.copy() * 0.0
        correction_odd = lastframe_odd.copy() * 0.0
        factor2_even = lastframe_even.copy() * 0.0
        factor2_odd = lastframe_odd.copy() * 0.0
        a1_even = lastframe_even.copy() * 0.0
        a1_odd = lastframe_odd.copy() * 0.0

        counts2_even = lastframe_even - crossopt_even
        counts2_odd = lastframe_odd - crossopt_odd

        counts2_even[np.where(counts2_even < 0)] = 0.0
        counts2_odd[np.where(counts2_odd < 0)] = 0.0

        # Find where counts2 > 0 and is finite
        good_even = np.where((counts2_even > 0) & np.isfinite(counts2_even))
        good_odd = np.where((counts2_odd > 0) & np.isfinite(counts2_odd))
        #__________________________________________________________________
        # even row values
        factor2_even[good_even] = 1.0 / \
            (np.exp(counts2_even[good_even] / b3_even) - 1)
        a1_even = b1_even * (np.power(counts2_even, b2_even)) * factor2_even
        #___________________________________________________________________
        # odd row values
        factor2_odd[good_odd] = 1.0 / \
            (np.exp(counts2_odd[good_odd] / b3_odd) - 1)
        a1_odd = b1_odd * (np.power(counts2_odd, b2_odd)) * factor2_odd
        #___________________________________________________________________
        # SATURATED DATA
        counts3_even = dn_lastfit[1::2, :] * sat_scale_even
        counts3_odd = dn_lastfit[0::2, :] * sat_scale_odd

        a1_sat_even = sat_final_slope_even * counts3_even + sat_mzp_even
        a1_sat_odd = sat_final_slope_odd * counts3_odd + sat_mzp_odd

        sat_even = sat[1::2, :]
        sat_odd = sat[0::2, :]

        # loop over groups in input science data:
        for j in range(sci_ngroups):

            # Compute the correction factors for even and odd rows
            T = (j + 1)
            eterm_even = np.exp(-T / tau_even)
            eterm_odd = np.exp(-T / tau_odd)

            # Apply the corrections to even and odd rows:
            # the first row is defined as odd (python index 0)
            # the second row is the first even row (python index of 1)
            correction_odd = lastframe_odd * a1_odd * 0.01 * eterm_odd
            correction_even = lastframe_even * a1_even * 0.01 * eterm_even
            correction_sat_odd = lastframe_odd * a1_sat_odd * 0.01 * eterm_odd
            correction_sat_even = lastframe_even * a1_sat_even * 0.01 * \
                eterm_even
            sat_index_even = np.where(sat_even)
            sat_index_odd = np.where(sat_odd)
            correction_even[sat_index_even] = \
                correction_sat_even[sat_index_even]
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

#    if subarray == 'SUBPRISM': subarray = 'SLITLESSPRISM'
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
    sat: the previous integration for this pixel saturated: yes/no
    dn_lastframe_23: extrapolated last frame using 2nd and 3rd to last frames
    dn_lastfrane_fit: extrapolated last frame using the fit to the entire ramp
    """

    nrows = input_model.data.shape[2]
    ncols = input_model.data.shape[3]
    dn_lastframe2 = input_model.data[i - 1][sci_ngroups - 2]
    dn_lastframe3 = input_model.data[i - 1][sci_ngroups - 3]
    dn_lastframe23 = dn_lastframe2.copy() * 0.0
    dn_lastframe_fit = dn_lastframe2.copy() * 0.0
    saturated = np.full((nrows, ncols), False)

    diff = dn_lastframe2 - dn_lastframe3
    dn_lastframe23 = dn_lastframe2 + diff

    # get saturation and reference pixel DQ flag values
    sat_flag = datamodels.dqflags.group['SATURATED']
    ref_flag = datamodels.dqflags.pixel['REFERENCE_PIXEL']

    # mark the locations of reference pixels
    refpix_2d = np.bitwise_and(input_model.pixeldq, ref_flag)
    dn_lastframe23[np.where(refpix_2d)] = 0.0

    # load the ramp data needed for computing slopes
    ramp3d = input_model.data[i - 1, 1:sci_ngroups - 1]
    groupdq3d = input_model.groupdq[i - 1, 1:sci_ngroups - 1]
    satmask3d = (groupdq3d == sat_flag)
    saturated = satmask3d.any(axis=0)

    # compute the slopes
    slope, intercept, ngood = ols_fit(ramp3d, groupdq3d)

    dn_lastframe_fit = slope * sci_ngroups + intercept

    # reset the results for pixels with zero slope
    slope0 = np.where(slope == 0)
    dn_lastframe_fit[slope0] = dn_lastframe23[slope0]

    # reset the results for reference pixels
    dn_lastframe23[np.where(refpix_2d)] = 0.0
    dn_lastframe_fit[np.where(refpix_2d)] = 0.0

    return saturated, dn_lastframe23, dn_lastframe_fit


def ols_fit(y, dq):

    """
    An estimation of the lastframe value from the previous integration is
    needed for the RSCD correction.
    This routine does a simple ordinary least squares fit to
    non-saturating data.
    """

    sat_flag = datamodels.dqflags.group['SATURATED']
    shape = y.shape

    # Find ramp values that are saturated
    x = np.arange(shape[0], dtype=np.float64)[:, np.newaxis, np.newaxis] * \
        np.ones(shape)
    good_data = np.bitwise_and(dq, sat_flag) == 0
    ngood = good_data.sum(axis=0)

    # Compute sums of unsaturated (good) x/y values
    sumx = (x * good_data).sum(axis=0)
    sumy = (y * good_data).sum(axis=0)
    sumxy = (x * y * good_data).sum(axis=0)
    sumxx = (x * x * good_data).sum(axis=0)
    nelem = good_data.sum(axis=0)

    # Compute the slopes and intercepts
    denom = nelem * sumxx - sumx * sumx
    with np.errstate(invalid='ignore'):  # ignore division warnings
        slope = (nelem * sumxy - sumx * sumy) / denom
        intercept = (sumxx * sumy - sumx * sumxy) / denom

    # Reset results to zero for pixels having < 3 unsaturated values
    bad = np.where(ngood < 3)
    slope[bad] = 0.0
    intercept[bad] = 0.0

    return (slope, intercept, ngood)
