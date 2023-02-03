#
#  Module for the RSCD correction for MIRI science data
#

import numpy as np
import logging

from stdatamodels.jwst.datamodels import dqflags

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
    input_model: ~jwst.datamodels.RampModel
        science data to be corrected

    rscd_model: ~jwst.datamodels.RSCDModel
        rscd reference data

    type: string
        type of algorithm ['baseline' or 'enhanced']

    Returns
    -------
    output_model: ~jwst.datamodels.RampModel
        RSCD-corrected science data

    """

    # Retrieve the reference parameters for this exposure type
    param = get_rscd_parameters(input_model, rscd_model)

    if not bool(param):  # empty dictionary
        log.warning('READPATT, SUBARRAY combination not found in ref file: RSCD correction will be skipped')
        input_model.meta.cal_step.rscd = 'SKIPPED'
        return input_model

    if type == 'baseline':
        group_skip = param['skip']
        output = correction_skip_groups(input_model, group_skip)
    else:
        # enhanced algorithm is not enabled yet (updated code and validation needed)
        log.warning('Enhanced algorithm not support yet: RSCD correction will be skipped')
        input_model.meta.cal_step.rscd = 'SKIPPED'
        return input_model
        # decay function algorithm update needed
        # output = correction_decay_function(input_model, param)

    return output


def correction_skip_groups(input_model, group_skip):
    """
    Short Summary
    -------------
    Set the initial groups in integration to DO_NOT_USE to skip groups
    affected by RSCD effect

    Parameters
    ----------
    input_model: ~jwst.datamodels.RampModel
        science data to be corrected

    group_skip: int
        number of groups to skip at the beginning of the ramp

    Returns
    -------
    output_model: ~jwst.datamodels.RampModel
        RSCD-corrected science data
    """

    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]       # number of integrations
    sci_ngroups = input_model.data.shape[1]     # number of groups

    log.debug("RSCD correction using: nints=%d, ngroups=%d" %
              (sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # If ngroups <= group_skip+3, skip the flagging
    # the +3 is to ensure there is a slope to be fit including the flagging for
    # the last frame correction
    if sci_ngroups <= (group_skip + 3):
        log.warning("Too few groups to apply RSCD correction")
        log.warning("RSCD step will be skipped")
        output.meta.cal_step.rscd = 'SKIPPED'
        return output

    # If ngroups > group_skip+3, set all of the GROUPDQ in the first group to 'DO_NOT_USE'
    output.groupdq[1:, 0:group_skip:, :] = \
        np.bitwise_or(output.groupdq[1:, 0:group_skip, :, :], dqflags.group['DO_NOT_USE'])
    log.debug(f"RSCD Sub: adding DO_NOT_USE to GROUPDQ for the first {group_skip} groups")
    output.meta.cal_step.rscd = 'COMPLETE'

    return output


def correction_decay_function(input_model, param):
    """
    Short Summary
    -------------
    Applies rscd correction to science arrays
    The last frame value from the previous integration is calculated two ways:
    1. using second and third to last frames to extrapolated to the last frame
    2. using the non saturating data, fit the data and extrapolate to last frame
    Because of the uncertainty of how well effects in th early part of the
    integration are corrected in the previous integration (reset anomaly, rscd
    effects, persistence) the lastframe determined from the second and third to
    last frames is considered a better estimate than that derived from a fit to
    the ramp.
    The last frame derived from fitting the non-saturating data is used in the
    correction if the previous integration saturated. This fit is extrapolated
    past saturation to estimate what the total number of electrons would have
    been.

    This correction has different correction parameters depending on whether
    the pixel is from an even row or odd row. The first row is define as an odd
    row. This even/odd row effect is likely a result of the reset electronics
    (MIRI resets in row pairs).

    Parameters
    ----------
    input_model: ~jwst.datamodels.RampModel
        science data to be corrected

    param: dict
        parameters of correction

    Returns
    -------
    output_model: ~jwst.datamodels.RampModel
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
        # __________________________________________________________________
        # even row values
        factor2_even[good_even] = 1.0 / \
            (np.exp(counts2_even[good_even] / b3_even) - 1)
        a1_even = b1_even * (np.power(counts2_even, b2_even)) * factor2_even
        # ___________________________________________________________________
        # odd row values
        factor2_odd[good_odd] = 1.0 / \
            (np.exp(counts2_odd[good_odd] / b3_odd) - 1)
        a1_odd = b1_odd * (np.power(counts2_odd, b2_odd)) * factor2_odd
        # ___________________________________________________________________
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
    Store the parameters in a param dictionary

    Parameters
    ----------
    input_model: ~jwst.datamodels.RampModel
        science data to be corrected

    rscd_model: ~jwst.datamodels.RSCDModel
        rscd reference data

    Returns
    -------
    param: dict
        dictionary of parameters

    """

    # Reference file parameters held in dictionary: param
    param = {}

    # read in the type of data from the input model (FAST,SLOW,FULL,SUBARRAY)
    readpatt = input_model.meta.exposure.readpatt
    subarray = input_model.meta.subarray.name

    # Check for old values of the MIRI LRS slitless subarray name
    # in the science data and change to the new
    if subarray.upper() == 'SUBPRISM':
        subarray = 'SLITLESSPRISM'

    # read table 1: containing the number of groups to skip
    for tabdata in rscd_model.rscd_group_skip_table:
        subarray_table = tabdata['subarray']
        readpatt_table = tabdata['readpatt']
        group_skip_table = tabdata['group_skip']
        if subarray_table == subarray and readpatt_table == readpatt:
            param['skip'] = group_skip_table
            break

    # read table 2: General RSCD enhanced parameters
    for tabdata in rscd_model.rscd_gen_table:
        readpatt_gen = tabdata['readpatt']
        subarray_gen = tabdata['subarray']
        lower_cutoff_gen = tabdata['lower_cutoff']
        alpha_even_gen = tabdata['alpha_even']
        alpha_odd_gen = tabdata['alpha_even']
        if subarray_gen == subarray and readpatt_gen == readpatt:
            param['gen'] = {}
            param['gen']['lower_cutoff'] = lower_cutoff_gen
            param['gen']['lower_alpha_odd'] = alpha_odd_gen
            param['gen']['lower_alpha_even'] = alpha_even_gen
            break

    # read table 3: Enhanced RSCD integration 1 parameters
    for tabdata in rscd_model.rscd_int1_table:
        readpatt_int1 = tabdata['readpatt']
        subarray_int1 = tabdata['subarray']
        rows_int1 = tabdata['rows']
        a0_int1 = tabdata['a0']
        a1_int1 = tabdata['a1']
        a2_int1 = tabdata['a2']
        a3_int1 = tabdata['a3']
        if subarray_int1 == subarray and readpatt_int1 == readpatt:
            param['int1'] = {}
            param['int1']['even'] = {}
            param['int1']['odd'] = {}
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
            break

    # read table 4: Enhanced RSCD integration 2 parameters
    for tabdata in rscd_model.rscd_int2_table:
        readpatt_int2 = tabdata['readpatt']
        subarray_int2 = tabdata['subarray']
        rows_int2 = tabdata['rows']
        a0_int2 = tabdata['b0']
        a1_int2 = tabdata['b1']
        a2_int2 = tabdata['b2']
        a3_int2 = tabdata['b3']
        if subarray_int2 == subarray and readpatt_int2 == readpatt:
            param['int2'] = {}
            param['int2']['even'] = {}
            param['int2']['odd'] = {}
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
            break

    # read table 5: Enhanced RSCD integration 3 parameters
    for tabdata in rscd_model.rscd_int3_table:
        readpatt_int3 = tabdata['readpatt']
        subarray_int3 = tabdata['subarray']
        rows_int3 = tabdata['rows']
        a0_int3 = tabdata['c0']
        a1_int3 = tabdata['c1']
        a2_int3 = tabdata['c2']
        a3_int3 = tabdata['c3']
        if subarray_int3 == subarray and readpatt_int3 == readpatt:
            param['int3'] = {}
            param['int3']['even'] = {}
            param['int3']['odd'] = {}
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
            break

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
    input_model: ~jwst.datamodels.RampModel
    i: integration #
    sci_ngroups: number of frames/integration

    return values
    -------------
    sat: the previous integration for this pixel saturated: yes/no
    dn_lastframe_23: extrapolated last frame using 2nd and 3rd to last frames
    dn_lastframe_fit: extrapolated last frame using the fit to the entire ramp
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
    sat_flag = dqflags.group['SATURATED']
    ref_flag = dqflags.pixel['REFERENCE_PIXEL']

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

    sat_flag = dqflags.group['SATURATED']
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
