#! /usr/bin/env python
#
#  ramp_fit.py - calculate weighted mean of slope, based on Massimo
#                Robberto's "On the Optimal Strategy to fit MULTIACCUM
#                ramps in the presence of cosmic rays."
#                (JWST-STScI-0001490,SM-12; 07/25/08).   The derivation
#                is a generalization for >1 cosmic rays, calculating
#                the slope and variance of the slope for each section
#                of the ramp (in between cosmic rays). The intervals are
#                determined from the input data quality arrays.
#
# Note:
# In this module, comments on the 'first group','second group', etc are
#    1-based, unless noted otherwise.

import time
import logging
import numpy as np

from .. import datamodels
from ..datamodels import dqflags
from ..lib import pipe_utils

from . import gls_fit           # used only if algorithm is "GLS"
from . import utils


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

BUFSIZE = 1024 * 30000  # 30Mb cache size for data section


def ramp_fit(model, buffsize, save_opt, readnoise_model, gain_model,
             algorithm, weighting):
    """
    Extended Summary
    ----------------
    Calculate the count rate for each pixel in all data cube sections and all
    integrations, equal to the slope for all sections (intervals between
    cosmic rays) of the pixel's ramp divided by the effective integration time.
    The weighting parameter must currently be set to 'optim', to use the optimal
    weighting (paper by Fixsen, ref. TBA) will be used in the fitting; this is
    currently the only supported weighting scheme.

    Parameters
    ----------
    model: data model
        input data model, assumed to be of type RampModel

    buffsize: int
        size of data section (buffer) in bytes

    save_opt: boolean
       calculate optional fitting results

    readnoise_model: instance of data Model
        readnoise for all pixels

    gain_model: instance of gain model
        gain for all pixels

    algorithm: string
        'OLS' specifies that ordinary least squares should be used;
        'GLS' specifies that generalized least squares should be used.

    weighting: string
        'optimal' specifies that optimal weighting should be used;
         currently the only weighting supported.

    Returns
    -------
    new_model: Data Model object
        DM object containing a rate image averaged over all integrations in
        the exposure

    int_model: Data Model object
        DM object containing rate images for each integration in the exposure

    opt_model: RampFitOutputModel object or None
        DM object containing optional OLS-specific ramp fitting data for the
        exposure

    gls_opt_model: GLS_RampFitModel object or None
        Object containing optional GLS-specific ramp fitting data for the
        exposure
    """
    if algorithm == "GLS":
        new_model, int_model, gls_opt_model = gls_ramp_fit(model,
                                buffsize, save_opt,
                                readnoise_model, gain_model)
        opt_model = None
    else:
        new_model, int_model, opt_model = \
               ols_ramp_fit(model, buffsize, save_opt, readnoise_model, \
               gain_model, weighting)
        gls_opt_model = None

    # Update data units in output models
    new_model.meta.bunit_data = 'DN/s'
    new_model.meta.bunit_err = 'DN/s'
    if int_model is not None:
        int_model.meta.bunit_data = 'DN/s'
        int_model.meta.bunit_err = 'DN/s'

    return new_model, int_model, opt_model, gls_opt_model


def ols_ramp_fit(model, buffsize, save_opt, readnoise_model, gain_model,
                 weighting):
    """
    Extended Summary
    ----------------
    Fit a ramp using ordinary least squares. Calculate the count rate for each
    pixel in all data cube sections and all integrations, equal to the weighted
    slope for all sections (intervals between cosmic rays) of the pixel's ramp
    divided by the effective integration time.

    Parameters
    ----------
    model: data model
        input data model, assumed to be of type RampModel

    buffsize: int
        size of data section (buffer) in bytes

    save_opt: boolean
        calculate optional fitting results

    readnoise_model: instance of data Model
        readnoise for all pixels

    gain_model: instance of gain model
        gain for all pixels

    weighting: string
        'optimal' specifies that optimal weighting should be used; currently
        the only weighting supported.

    Returns
    -------
    new_model: Data Model object
        DM object containing a rate image averaged over all integrations in
        the exposure

    int_model: Data Model object or None
        DM object containing rate images for each integration in the exposure,
        or None if there is only one integration in the exposure

    opt_model: Data Model object or None
        DM object containing optional OLS-specific ramp fitting data for the
        exposure; this will be None if save_opt is False
    """
    tstart = time.time()

    # Get needed sizes and shapes
    nreads, npix, imshape, cubeshape, n_int, instrume, frame_time, ngroups, \
        group_time = utils.get_dataset_info(model)

    # Save original shapes for writing to log file, as these may change for MIRI
    orig_nreads = nreads
    orig_cubeshape = cubeshape

    # For MIRI datasets having >1 groups, if all final groups are flagged as 
    # DO_NOT_USE, resize the input model arrays to exclude the final group.
    # Similarly, if all first groups are flagged as DO_NOT_USE, resize the input 
    # model arrays to exclude the first group.

    if (instrume == 'MIRI' and nreads > 1):
        first_gdq = model.groupdq[:,1,:,:]

        if np.all(np.bitwise_and( first_gdq, dqflags.group['DO_NOT_USE'] )):
            model.data = model.data[:,1:,:,:]
            model.err = model.err[:,1:,:,:]
            model.groupdq = model.groupdq[:,1:,:,:]
            nreads -= 1
            ngroups -= 1
            cubeshape = (nreads,)+imshape
            log.info('MIRI dataset has all first groups flagged as DO_NOT_USE.')

        last_gdq = model.groupdq[:,-1,:,:]
        if np.all(np.bitwise_and( last_gdq, dqflags.group['DO_NOT_USE'] )):
            model.data = model.data[:,:-1,:,:]
            model.err = model.err[:,:-1,:,:]
            model.groupdq = model.groupdq[:,:-1,:,:]
            nreads -= 1
            ngroups -= 1
            cubeshape = (nreads,)+imshape
            log.info('MIRI dataset has all final groups flagged as DO_NOT_USE.')

        # Next block is to satisfy github issue 1681:
        # "MIRI FirstFrame and LastFrame minimum number of groups #1681"
        if (ngroups < 2):
            log.warning('MIRI datasets require at least 2 groups/integration')
            log.warning('(NGROUPS), so will not process this dataset.')
            return None, None, None

    if (ngroups == 1):
        log.warning('Dataset has NGROUPS=1, so count rates for each integration')
        log.warning('will be calculated as the value of that 1 group divided by')
        log.warning('the group exposure time.')

    # Calculate effective integration time (once EFFINTIM has been populated
    #   and accessible, will use that instead), and other keywords that will
    #   needed if the pedestal calculation is requested. Note 'nframes'
    #   is the number of given by the NFRAMES keyword, and is the number of
    #   frames averaged on-board for a group, i.e., it does not include the
    #   groupgap.
    effintim, nframes, groupgap, dropframes1= utils.get_efftim_ped(model)

    # Get GROUP DQ and ERR arrays from input file
    gdq_cube = model.groupdq
    gdq_cube_shape = gdq_cube.shape

    # Get max number of segments fit in all integrations
    max_seg = calc_num_seg(gdq_cube, n_int)
    del gdq_cube

    f_max_seg = 0  # final number to use, usually overwritten by actual value

    (dq_int, median_diffs_2d, num_seg_per_int, sat_0th_group_int) =\
        utils.alloc_arrays_1(n_int, imshape)

    opt_res = utils.OptRes(n_int, imshape, max_seg, nreads, save_opt)

    # Get Pixel DQ array from input file. The incoming RampModel has uint8
    #   PIXELDQ, but ramp fitting will update this array here by flagging
    #   the 2D PIXELDQ locations where the ramp data has been previously
    #   flagged as jump-detected or saturated. These additional bit values
    #   require this local variable to be uint16, and it will be used as the
    #   (uint16) PIXELDQ in the outgoing ImageModel.
    pixeldq = model.pixeldq.copy()

    # calculate number of (contiguous) rows per data section
    nrows = calc_nrows(model, buffsize, cubeshape, nreads)

    # Get readnoise array for calculation of variance of noiseless ramps, and
    #   gain array in case optimal weighting is to be done
    readnoise_2d, gain_2d = utils.get_ref_subs(model, readnoise_model,
                                               gain_model, nframes)

    # Flag any bad pixels in the gain
    pixeldq = utils.reset_bad_gain( pixeldq, gain_2d )

    # In this 'First Pass' over the data, loop over integrations and data
    #   sections to calculate the estimated median slopes, which will be used
    #   to calculate the variances. This is the same method to estimate slopes
    #   as is done in the jump detection step, except here CR-affected and
    #   saturated groups have already been flagged. The actual, fit, slopes for
    #   each segment are also calculated here.

    # Loop over data integrations:
    for num_int in range(0, n_int):
        # Loop over data sections
        for rlo in range(0, cubeshape[1], nrows):
            rhi = rlo + nrows

            if rhi > cubeshape[1]:
                rhi = cubeshape[1]

            data_sect = model.get_section('data')[num_int, :, rlo:rhi, :]

            # Skip data section if it is all NaNs
            if  np.all(np.isnan( data_sect)):
                log.error('Current data section is all nans, so not processing the section.')
                continue

            # first frame section for 1st group of current integration
            ff_sect = model.get_section('data')[ num_int, 0, rlo:rhi, :].\
                astype(np.float32)

            # Get appropriate sections
            gdq_sect = model.get_section('groupdq')[num_int, :, rlo:rhi, :]

            rn_sect = readnoise_2d[rlo:rhi, :]
            gain_sect = gain_2d[rlo:rhi, :]

            # Reset all saturated groups in the input data array to NaN
            where_sat = np.where( np.bitwise_and(gdq_sect,
                                  dqflags.group['SATURATED']) != 0)

            data_sect[ where_sat ] = np.NaN
            del where_sat

            # Compute the first differences of all groups
            first_diffs_sect = np.diff(data_sect, axis=0)

            # If the dataset has only 1 group/integ, assume the 'previous group'
            #   is all zeros, so just use data as the difference
            if (first_diffs_sect.shape[0] == 0):
                first_diffs_sect = data_sect.copy()
            else:
                # Similarly, for datasets having >1 group/integ and having
                #  single-group segments, just use the data as the difference
                wh_nan = np.where( np.isnan( first_diffs_sect[0,:,:]) )

                if (len(wh_nan[0]) > 0):
                    first_diffs_sect[0,:,:][wh_nan] = data_sect[0,:,:][wh_nan]

                del wh_nan

                i_group,i_yy,i_xx, = np.where( np.bitwise_and(gdq_sect[1:,:,:],
                                              dqflags.group['JUMP_DET']) != 0)

                # Mask all the first differences that are affected by a CR,
                #   starting at group 1.  The purpose of starting at index 1 is
                #   to shift all the indices down by 1, so they line up with the
                #   indices in first_diffs.
                first_diffs_sect[ i_group-1, i_yy, i_xx ] = np.NaN

                del i_group, i_yy, i_xx

                # Check for pixels in which there is good data in 0th group, but
                #   all first_diffs for this ramp are NaN because there are too
                #   few good groups past the 0th. Due to the shortage of good
                #   data, the first_diffs will be set here equal to the data in
                #   the 0th group.
                wh_min = np.where( np.logical_and(np.isnan(first_diffs_sect)
                              .all(axis=0), np.isfinite( data_sect[0,:,:])))
                if len(wh_min[0] > 0):
                    first_diffs_sect[0,:,:][wh_min] = data_sect[0,:,:][wh_min]

                del wh_min

            # All first differences affected by saturation and CRs have been set
            #  to NaN, so compute the median of all non-NaN first differences.
            nan_med = np.nanmedian(first_diffs_sect, axis=0)
            nan_med[np.isnan(nan_med)] = 0. # if all first_diffs_sect are nans
            median_diffs_2d[ rlo:rhi, : ] += nan_med

            # Calculate the slope of each segment
            t_dq_cube, inv_var, opt_res, f_max_seg, num_seg = \
                 calc_slope(data_sect, gdq_sect, frame_time, opt_res, save_opt,
                            rn_sect, gain_sect, max_seg, ngroups, weighting,
                            f_max_seg, group_time)
           
            del gain_sect

            # Populate 3D num_seg { integ, y, x } with 2D num_seg for this data
            #  section (y,x) and integration (num_int)
            sect_shape = data_sect.shape[-2:]
            num_seg_per_int[num_int, rlo:rhi, :] = num_seg.reshape(sect_shape)

            # Populate integ-spec slice which is set if 0th group has SAT
            wh_sat0 = np.where( np.bitwise_and(gdq_sect[0,:,:], 
                                dqflags.group['SATURATED']))
            if ( len(wh_sat0[0] ) > 0):
                sat_0th_group_int[num_int, rlo:rhi, :][ wh_sat0 ] = 1

            del wh_sat0

            pixeldq_sect = pixeldq[rlo:rhi, :].copy()
            dq_int[num_int, rlo:rhi, :] = \
                  dq_compress_sect(t_dq_cube, pixeldq_sect).copy()

            del t_dq_cube

            # Loop over the segments and copy the reshaped 2D segment-specific
            #  results for the current data section to the 4D output arrays.
            opt_res.reshape_res(num_int, rlo, rhi, sect_shape, ff_sect, save_opt)

            if save_opt:
                # Calculate difference between each slice and the previous slice
                #   as approximation to cosmic ray amplitude for those pixels
                #   having their DQ set for cosmic rays
                data_diff = data_sect - utils.shift_z(data_sect, -1)
                dq_cr = np.bitwise_and(dqflags.group['JUMP_DET'], gdq_sect)

                opt_res.cr_mag_seg[num_int, :, rlo:rhi, :] = \
                       data_diff * (dq_cr != 0)

                del data_diff
                
            del data_sect
            del ff_sect
            del gdq_sect
  
    if pixeldq_sect is not None:
        del pixeldq_sect

    # Compute the final 2D array of differences; create rate array
    median_diffs_2d /= n_int
    med_rates = median_diffs_2d/group_time

    del median_diffs_2d
    del first_diffs_sect

    (var_p3, var_r3, var_p4, var_r4, var_both4, var_both3,
     inv_var_both4, s_inv_var_p3, s_inv_var_r3, s_inv_var_both3,
     segs_4) = utils.alloc_arrays_2(n_int, imshape, max_seg)

    # In this 'Second Pass' over the data, loop over integrations and data
    #   sections to calculate the variances of the slope using the estimated
    #   median slopes from the 'First Pass'. These variances are due to Poisson
    #   noise only, read noise only, and the combination of Poisson noise and
    #   read noise. The integration-specific variances are 3D arrays, and
    #   the segment-specific variances are 4D arrays . The naming convention for
    #   the arrays:
    #     'var': a variance
    #     'p3': intermediate 3D array for variance due to Poisson noise
    #     'r4': intermediate 4D array for variance due to read noise
    #     'both4': intermediate 4D array for combined variance due to both
    #               Poisson and read noise
    #     'inv_<X>': intermediate array = 1/<X>
    #     's_inv_<X>': intermediate array = 1/<X>, summed over integrations

    # Loop over data integrations
    for num_int in range(n_int):

        # Loop over data sections
        for rlo in range(0, cubeshape[1], nrows):
            rhi = rlo + nrows

            if rhi > cubeshape[1]:
                rhi = cubeshape[1]

            # data_sect and gdq_sect are: [ groups, y, x ]
            data_sect = model.get_section('data')[num_int, :, rlo:rhi, :]
            gdq_sect = model.get_section('groupdq')[num_int, :, rlo:rhi, :]

            rn_sect = readnoise_2d[rlo:rhi, :]
            gain_sect = gain_2d[rlo:rhi, :]

            # Calculate results needed to compute the variance arrays
            den_r3, den_p3, num_r3, segs_beg_3 = \
                utils.calc_slope_vars( rn_sect, gain_sect, gdq_sect,
                                       group_time, max_seg )

            segs_4[num_int, :, rlo:rhi, :] = segs_beg_3
            var_p4[num_int, :, rlo:rhi, :] = den_p3 * med_rates[ rlo:rhi, :]
            var_r4[num_int, :, rlo:rhi, :] = num_r3 * den_r3

            del den_r3, den_p3, num_r3, segs_beg_3
            del gain_sect
            del gdq_sect

        # The next 4 statements zero out entries for non-existing segments, and
        #   set the variances for segments having negative slopes (the segment
        #   variance is proportional to the median estimated slope) to
        #   outrageously large values so that they will have negligible
        #   contributions.
        var_p4[num_int,:,:,:] *= ( segs_4[num_int,:,:,:] > 0)
        var_p4[var_p4 <= 0.] = utils.LARGE_VARIANCE
        var_r4[num_int,:,:,:] *= ( segs_4[num_int,:,:,:] > 0)
        var_r4[var_r4 <= 0.] = utils.LARGE_VARIANCE

        s_inv_var_p3[num_int, :, :] = (1./var_p4[num_int, :, :, :]).sum(axis=0)
        var_p3[num_int, :, :] = 1./ s_inv_var_p3[num_int, :, :]
        s_inv_var_r3[num_int, :, :] = (1./var_r4[num_int, :, :, :]).sum(axis=0)
        var_r3[num_int, :, :] = 1./ s_inv_var_r3[num_int, :, :]

        var_both4[ num_int,:,:,:] = var_r4[num_int,:,:,:] + var_p4[num_int,:,:,:]
        inv_var_both4[num_int, :, :, :] = 1./var_both4[num_int, :, :, :]

        # Want to retain values in the 4D arrays only for the segments that each
        #   pixel has, so will zero out values for the higher indices. Creating
        #   and manipulating intermediate arrays (views, such as var_p4_int
        #   will zero out the appropriate indices in var_p4 and var_r4.)
        # Extract the slice of 4D arrays for the current integration
        #
        var_p4_int = var_p4[ num_int,:,:,:]   # [ segment, y, x ]
        inv_var_both4_int = inv_var_both4[ num_int,:,:,:]

        # Zero out non-existing segments
        var_p4_int *= ( segs_4[num_int,:,:,:] > 0)
        inv_var_both4_int *= ( segs_4[num_int,:,:,:] > 0)

        # reshape these arrays to simplify masking [ segment, 1D pixel ]
        var_p4_int2 = var_p4_int.reshape(( var_p4_int.shape[0],
                            var_p4_int.shape[1]*var_p4_int.shape[2]))

        num_seg_int = num_seg_per_int[num_int,:, :] # number of segs per pixel

        s_inv_var_both3[num_int,:,:] = (inv_var_both4[num_int,:,:,:]).sum(axis=0)
        var_both3[num_int, :, :] = 1./s_inv_var_both3[num_int, :, :]

        del var_p4_int
        del var_p4_int2

    del gain_2d

    var_p4 *= ( segs_4[:,:,:,:] > 0) # Zero out non-existing segments
    var_r4 *= ( segs_4[:,:,:,:] > 0)

    # Delete lots of arrays no longer needed
    if inv_var_both4_int is not None:
        del inv_var_both4_int

    if med_rates is not None:
        del med_rates

    if num_seg_per_int is not None:
        del num_seg_per_int

    if readnoise_2d is not None:
        del readnoise_2d

    if rn_sect is not None:
        del rn_sect

    if segs_4 is not None:
        del segs_4

    # Now that the segment-specific and integration-specific variances have
    #   been calculated, the segment-specific, integration-specific, and
    #   overall slopes will be calculated. The integration-specific slope is
    #   calculated as a weighted average of the segments in the integration:
    #     slope_int = sum_over_segs(slope_seg/var_seg)/ sum_over_segs(1/var_seg)
    #  The overall slope is calculated as a weighted average of the segments in
    #     all integrations:
    #     slope = sum_over_integs_and_segs(slope_seg/var_seg)/
    #                    sum_over_integs_and_segs(1/var_seg)
    slope_by_var4 = opt_res.slope_seg.copy()/var_both4 

    del var_both4

    s_slope_by_var3 = slope_by_var4.sum(axis=1) # sum over segments (not integs)
    s_slope_by_var2 = s_slope_by_var3.sum(axis=0) # sum over integrations
    s_inv_var_both2 = s_inv_var_both3.sum(axis=0)

    # Compute the 'dataset-averaged' slope
    slope_dataset2 = s_slope_by_var2/s_inv_var_both2

    del s_inv_var_both2, s_slope_by_var2, s_slope_by_var3, slope_by_var4
    del s_inv_var_both3

    #  Replace nans in slope_dataset2 with 0 (for non-existing segments)
    slope_dataset2[np.isnan(slope_dataset2)] = 0.

    # Compute the integration-specific slope
    the_num = (opt_res.slope_seg * inv_var_both4).sum(axis = 1)

    the_den = (inv_var_both4).sum(axis=1)
    slope_int = the_num/the_den

    del the_num, the_den

    # Clean up ramps that are SAT on their initial groups; set ramp parameters
    #   for variances and slope so they will not contribute
    var_p3, var_both3, slope_int = utils.fix_sat_ramps( sat_0th_group_int,
                                           var_p3, var_both3, slope_int)

    if sat_0th_group_int is not None:
        del sat_0th_group_int

    # Loop over data integrations to calculate integration-specific pedestal
    if save_opt:
        dq_slice = np.zeros((gdq_cube_shape[2],gdq_cube_shape[3]),
                            dtype=np.uint32)

        for num_int in range(0, n_int):
            dq_slice =  model.get_section('groupdq')[num_int, 0, :, :]
            opt_res.ped_int[ num_int, :, : ] = \
                utils.calc_pedestal(num_int, slope_int, opt_res.firstf_int,
                    dq_slice, nframes, groupgap, dropframes1)

        del dq_slice

    # Collect optional results for output
    if save_opt:
        gdq_cube = model.groupdq  
        opt_res.shrink_crmag(n_int, gdq_cube, imshape, nreads)
        del gdq_cube 

        # Truncate results at the maximum number of segments found
        opt_res.slope_seg = opt_res.slope_seg[:,:f_max_seg,:,:]
        opt_res.sigslope_seg = opt_res.sigslope_seg[:,:f_max_seg,:,:]
        opt_res.yint_seg = opt_res.yint_seg[:,:f_max_seg,:,:]
        opt_res.sigyint_seg = opt_res.sigyint_seg[:,:f_max_seg,:,:]
        opt_res.weights = (inv_var_both4[:,:f_max_seg,:,:])**2.
        opt_res.var_p_seg = var_p4[:,:f_max_seg,:,:]
        opt_res.var_r_seg = var_r4[:,:f_max_seg,:,:]

        opt_model = opt_res.output_optional(model, effintim)
    else:
        opt_model = None

    if inv_var_both4 is not None:
        del inv_var_both4

    if var_p4 is not None:
        del var_p4

    if var_r4 is not None:
        del var_r4

    if inv_var is not None:
        del inv_var

    if pixeldq is not None:  
        del pixeldq

    # For multiple-integration datasets, will output integration-specific
    #    results to separate file named <basename> + '_integ.fits'
    int_times = None
    if n_int > 1:
        if pipe_utils.is_tso(model) and hasattr(model, 'int_times'):
            int_times = model.int_times
        else:
            int_times = None
        int_model = utils.output_integ(model, slope_int, dq_int, effintim,
                                       var_p3, var_r3, var_both3, int_times)
    else:
        int_model = None

    if opt_res is not None:
        del opt_res

    if slope_int is not None:
        del slope_int
    del var_p3
    del var_r3
    del var_both3
    if int_times is not None:
        del int_times

    # Divide slopes by total (summed over all integrations) effective
    #   integration time to give count rates.
    c_rates = slope_dataset2 / effintim

    # Compress all integration's dq arrays to create 2D PIXELDDQ array for
    #   primary output
    final_pixeldq = dq_compress_final(dq_int, n_int)

    if dq_int is not None:
        del dq_int

    tstop = time.time()

    log_stats(c_rates)

    log.debug('Instrument: %s', instrume)
    log.debug('Number of pixels in 2D array: %d', npix)
    log.debug('Shape of 2D image: (%d, %d)' %(imshape))
    log.debug('Shape of data cube: (%d, %d, %d)' %(orig_cubeshape))
    log.debug('Buffer size (bytes): %d', buffsize)
    log.debug('Number of rows per buffer: %d', nrows)
    log.info('Number of groups per integration: %d', orig_nreads)
    log.info('Number of integrations: %d', n_int)
    log.debug('The execution time in seconds: %f', tstop - tstart)

    # Compute 2D 'error' for primary output; this is the standard deviation due
    # to both the Poisson and read noise
    var_p2 = 1/(s_inv_var_p3.sum(axis=0))
    var_r2 = 1/(s_inv_var_r3.sum(axis=0))

    del s_inv_var_p3
    del s_inv_var_r3

    # Create new model...
    new_model = datamodels.ImageModel(data=c_rates.astype(np.float32),
            dq=final_pixeldq.astype(np.uint32),
            var_poisson=var_p2.astype(np.float32),
            var_rnoise=var_r2.astype(np.float32), err=np.sqrt(var_p2 + var_r2))

    new_model.update(model)  # ... and add all keys from input

    return new_model, int_model, opt_model


def gls_ramp_fit(model,
                 buffsize, save_opt,
                 readnoise_model, gain_model):
    """Fit a ramp using generalized least squares.

    Extended Summary
    ----------------
    Calculate the count rate for each pixel in the data ramp, for every
    integration.  Generalized least squares is used for fitting the ramp
    in order to take into account the correlation between reads.  If the
    input file contains multiple integrations, a second output file will
    be written, containing per-integration count rates.

    One additional file can optionally be written (if save_opt is True),
    containing per-integration data.

    Parameters
    ----------
    model: data model
        Input data model, assumed to be of type RampModel.

    buffsize: int
        Size of data section (buffer) in bytes.

    save_opt: boolean
        Calculate optional fitting results.

    readnoise_model: instance of data Model
        Readnoise for all pixels.

    gain_model: instance of gain model
        Gain for all pixels.

    Returns
    -------
    new_model: Data Model object
        DM object containing a rate image averaged over all integrations in
        the exposure.

    int_model: Data Model object or None
        DM object containing rate images for each integration in the exposure,
        or None if there is only one integration.

    gls_opt_model: GLS_RampFitModel object or None
        Object containing optional GLS-specific ramp fitting data for the
        exposure; this will be None if save_opt is False.
    """

    tstart = time.time()

    # get needed sizes and shapes
    nreads, npix, imshape, cubeshape, n_int, instrume, frame_time, ngroups = \
            utils.get_dataset_info(model)

    (group_time, nframes_used, saturated_flag, jump_flag) = \
            utils.get_more_info(model)
    if n_int > 1:
        # `slopes` will be used for accumulating the sum of weighted slopes.
        med_slopes = np.zeros(imshape, dtype=np.float64)
        sum_weight = np.zeros(imshape, dtype=np.float64)

    # For multiple-integration datasets, will output integration-specific
    # results to separate file named <basename> + '_rateints.fits'.
    # Even if there's only one integration, the output results will be
    # saved in these arrays.
    slope_int = np.zeros((n_int,) + imshape, dtype=np.float32)
    slope_err_int = np.zeros((n_int,) + imshape, dtype=np.float32)
    dq_int = np.zeros((n_int,) + imshape, dtype=np.uint32)

    # Get GROUP DQ array from input file
    gdq_cube = model.groupdq

    # Determine the maximum number of cosmic ray hits for any pixel.
    max_num_cr = -1                     # invalid initial value
    for num_int in range(n_int):
        i_max_num_cr = utils.get_max_num_cr(gdq_cube[num_int,:,:,:], jump_flag)
        max_num_cr = max(max_num_cr, i_max_num_cr)

    if save_opt:
        # Create arrays for the fitted values of zero-point intercept and
        # cosmic-ray amplitudes, and their errors.
        intercept_int = np.zeros((n_int,) + imshape, dtype=np.float32)
        intercept_err_int = np.zeros((n_int,) + imshape, dtype=np.float32)
        # The pedestal is the extrapolation of the first group back to zero
        # time, for each integration.
        pedestal_int = np.zeros((n_int,) + imshape, dtype=np.float32)
        # The first group, for calculating the pedestal.  (This only needs
        # to be nrows high, but we don't have nrows yet.  xxx)
        first_group = np.zeros(imshape, dtype=np.float32)
        # If there are no cosmic rays, set the last axis length to 1.
        shape_ampl = (n_int, imshape[0], imshape[1], max(1, max_num_cr))
        ampl_int = np.zeros(shape_ampl, dtype=np.float32)
        ampl_err_int = np.zeros(shape_ampl, dtype=np.float32)

    # Used for flagging pixels with UNRELIABLE_SLOPE.
    temp_dq = np.zeros(imshape, dtype=np.uint32)

    # Get Pixel DQ array from input file. The incoming RampModel has uint8
    # PIXELDQ, but ramp fitting will update this array here by flagging
    # the 2D PIXELDQ locations where the ramp data has been previously
    # flagged as jump-detected or saturated. These additional bit values
    # require this local variable to be uint16, and it will be used as the
    # (uint16) PIXELDQ in the outgoing ImageModel.
    pixeldq = model.pixeldq.copy()

    # calculate number of (contiguous) rows per data section
    nrows = calc_nrows(model, buffsize, cubeshape, nreads)

    # Get readnoise array and gain array
    readnoise_2d, gain_2d = utils.get_ref_subs(model, readnoise_model,
                                               gain_model)

    # Flag any bad pixels in the gain
    pixeldq = utils.reset_bad_gain( pixeldq, gain_2d )

    # loop over data integrations
    for num_int in range(n_int):
        if save_opt:
            first_group[:, :] = 0.      # re-use this for each integration

        # loop over data sections
        for rlo in range(0, cubeshape[1], nrows):
            rhi = rlo + nrows

            if rhi > cubeshape[1]:
                rhi = cubeshape[1]

            data_sect = model.get_section('data')[num_int, :, rlo:rhi, :]

            # We'll propagate error estimates from previous steps to the
            # current step by using the variance.
            input_var_sect = model.get_section('err')[num_int, :, rlo:rhi, :]
            input_var_sect = input_var_sect**2

            gdq_sect = gdq_cube[num_int, :, rlo:rhi, :]
            rn_sect = readnoise_2d[rlo:rhi, :]
            gain_sect = gain_2d[rlo:rhi, :]

            # Convert the data section from DN to electrons.
            data_sect *= gain_sect
            if save_opt:
                first_group[rlo:rhi, :] = data_sect[0, :, :].copy()

            (intercept_sect, intercept_var_sect,
             slope_sect, slope_var_sect,
             cr_sect, cr_var_sect) = \
                        gls_fit.determine_slope(data_sect, input_var_sect,
                                                gdq_sect, rn_sect, gain_sect,
                                                frame_time, group_time,
                                                nframes_used, max_num_cr,
                                                saturated_flag, jump_flag)

            slope_int[num_int, rlo:rhi, :] = slope_sect.copy()
            v_mask = (slope_var_sect <= 0.)
            if v_mask.any():
                # Replace negative or zero variances with a large value.
                slope_var_sect[v_mask] = LARGE_VARIANCE
                # Also set a flag in the pixel dq array.
                temp_dq[rlo:rhi, :][v_mask] = dqflags.pixel['UNRELIABLE_SLOPE']
            del v_mask
            # If a pixel was flagged (by an earlier step) as saturated in
            # the first group, flag the pixel as bad.
            # Note:  save s_mask until after the call to utils.gls_pedestal.
            s_mask = (gdq_sect[0] == saturated_flag)
            if s_mask.any():
                temp_dq[rlo:rhi, :][s_mask] = dqflags.pixel['UNRELIABLE_SLOPE']
            slope_err_int[num_int, rlo:rhi, :] = np.sqrt(slope_var_sect)

            # We need to take a weighted average if (and only if) n_int > 1.
            # Accumulate sum of slopes and sum of weights.
            if n_int > 1:
                weight = 1. / slope_var_sect
                slopes[rlo:rhi, :] += (slope_sect * weight)
                sum_weight[rlo:rhi, :] += weight

            if save_opt:
                # Save the intercepts and cosmic-ray amplitudes for the
                # current integration.
                intercept_int[num_int, rlo:rhi, :] = intercept_sect.copy()
                intercept_err_int[num_int, rlo:rhi, :] = \
                        np.sqrt(np.abs(intercept_var_sect))
                pedestal_int[num_int, rlo:rhi, :] = \
                        utils.gls_pedestal(first_group[rlo:rhi, :],
                                           slope_int[num_int, rlo:rhi, :],
                                           s_mask,
                                           frame_time, nframes_used)
                ampl_int[num_int, rlo:rhi, :, :] = cr_sect.copy()
                ampl_err_int[num_int, rlo:rhi, :, :] = \
                        np.sqrt(np.abs(cr_var_sect))
            del s_mask

            # Compress 4D->2D dq arrays for saturated and jump-detected
            #   pixels
            pixeldq_sect = pixeldq[rlo:rhi, :].copy()
            dq_int[num_int, rlo:rhi, :] = \
                  dq_compress_sect(gdq_sect, pixeldq_sect).copy()

        # temp_dq |= dq_int[num_int, :, :]
        # dq_int[num_int, :, :] = temp_dq.copy()
        dq_int[num_int, :, :] |= temp_dq
        temp_dq[:, :] = 0               # initialize for next integration

    # Average the slope over all integrations.
    if n_int > 1:
        sum_weight = np.where(sum_weight <= 0., 1., sum_weight)
        recip_sum_weight = 1. / sum_weight
        slopes *= recip_sum_weight
        gls_err = np.sqrt(recip_sum_weight)

    # Convert back from electrons to DN.
    slope_int /= gain_2d
    slope_err_int /= gain_2d
    if n_int > 1:
        slopes /= gain_2d
        gls_err /= gain_2d
    if save_opt:
        intercept_int /= gain_2d
        intercept_err_int /= gain_2d
        pedestal_int /= gain_2d
        gain_shape = gain_2d.shape
        gain_4d = gain_2d.reshape((1, gain_shape[0], gain_shape[1], 1))
        ampl_int /= gain_4d
        ampl_err_int /= gain_4d
        del gain_4d
    del gain_2d

    # Compress all integration's dq arrays to create 2D PIXELDDQ array for
    #   primary output
    final_pixeldq = dq_compress_final(dq_int, n_int)

    if n_int > 1:
        effintim = 1.                   # slopes are already in DN/s
        int_model = utils.output_integ(model, slope_int, slope_err_int, dq_int,
                                       effintim)
    else:
        int_model = None

    if save_opt: # collect optional results for output
        # Get the zero-point intercepts and the cosmic-ray amplitudes for
        # each integration (even if there's only one integration).
        gls_opt_model = utils.gls_output_optional(model,
                                intercept_int, intercept_err_int,
                                pedestal_int,
                                ampl_int, ampl_err_int)
    else:
        gls_opt_model = None

    tstop = time.time()

    if n_int > 1:
        log_stats(slopes)
    else:
        log_stats(slope_int[0])

    log.info('Instrument: %s' % instrume)
    log.info('Number of groups per integration: %d' % nreads)
    log.info('Number of pixels in 2D array: %d' % npix)
    log.info('Shape of 2D image: (%d, %d)' % imshape)
    log.info('Shape of data cube: (%d, %d, %d)' % cubeshape)
    log.info('Buffer size (bytes): %d' % buffsize)
    log.info('Number of rows per buffer: %d' % nrows)
    log.info('Number of integrations: %d' % n_int)
    log.info('The execution time in seconds: %f' % (tstop - tstart,))

    # Create new model...
    if n_int > 1:
        new_model = datamodels.ImageModel(data=slopes.astype(np.float32),
                           dq=final_pixeldq, err=gls_err.astype(np.float32))
    else:
        new_model = datamodels.ImageModel(data=slope_int[0], dq=final_pixeldq,
                                          err=slope_err_int[0])

    new_model.update(model)     # ... and add all keys from input

    return new_model, int_model, gls_opt_model


def calc_power(snr):
    """
    Short Summary
    -------------
    Using the given SNR, calculate the weighting exponent, which is from
    `Fixsen, D.J., Offenberg, J.D., Hanisch, R.J., Mather, J.C, Nieto,
    Santisteban, M.A., Sengupta, R., & Stockman, H.S., 2000, PASP, 112, 1350`.

    Parameters
    ----------
    snr: float, 1D
        signal-to-noise for the ramp segments

    Returns
    -------
    pow_wt.ravel(): float, 1D
        weighting exponent
    """
    pow_wt = snr.copy() * 0.0
    pow_wt[np.where(snr > 5.)] = 0.4
    pow_wt[np.where(snr > 10.)] = 1.0
    pow_wt[np.where(snr > 20.)] = 3.0
    pow_wt[np.where(snr > 50.)] = 6.0
    pow_wt[np.where(snr > 100.)] = 10.0

    return pow_wt.ravel()


def dq_compress_final(dq_int, n_int):
    """
    Extended Summary
    ----------------
    Combine the integration-specific dq arrays (which have already been
    compressed and combined with the PIXELDQ array) to create the dq array
    of the primary output product.

    Parameters
    ----------
    dq_int: uint16, 3D array
        cube of combined dq arrays for all data sections in a single itegration

    n_int: int
        total number of integrations in data set

    Returns
    -------
    f_dq: float, 2D array
        combination of all integration's pixeldq arrays

    """
    f_dq = dq_int[0, :, :]

    for jj in range(1, n_int):
        f_dq = np.bitwise_or(f_dq, dq_int[jj, :, :])

    return f_dq


def dq_compress_sect(gdq_sect, pixeldq_sect):
    """
    Extended Summary
    ----------------
    Get ramp locations where the data has been flagged as saturated in the 4D
    GROUPDQ array for the current data section, find the corresponding image
    locations, and set the SATURATED flag in those locations in the PIXELDQ
    array. Similarly, get the ramp locations where the data has been flagged as
    a jump detection in the 4D GROUPDQ array, find the corresponding image
    locations, and set the COSMIC_BEFORE flag in those locations in the PIXELDQ
    array. These modifications to the section of the PIXELDQ array are not used
    to flag groups for any computations; they are used only in the integration-
    specific output.

    Parameters
    ----------
    gdq_sect: int (uint8), 3D array
        cube of GROUPDQ array for a data section

    pixeldq_sect: int (uint16), 2D array
        dq array of data section of input model

    Returns
    -------
    pixeldq_sect: int (uint16), 2D array
        dq array of data section updated with saturated and jump-detected flags

    """
    sat_loc_r = np.bitwise_and(gdq_sect, dqflags.group['SATURATED'])
    sat_loc_im = np.where(sat_loc_r.sum(axis=0) > 0)
    pixeldq_sect[sat_loc_im] = np.bitwise_or(pixeldq_sect[sat_loc_im],
                                             dqflags.pixel['SATURATED'])

    cr_loc_r = np.bitwise_and(gdq_sect, dqflags.group['JUMP_DET'])
    cr_loc_im = np.where(cr_loc_r.sum(axis=0) > 0)
    pixeldq_sect[cr_loc_im] = np.bitwise_or(pixeldq_sect[cr_loc_im],
                                            dqflags.pixel['JUMP_DET'])

    return pixeldq_sect


def calc_nrows(model, buffsize, cubeshape, nreads):
    """
    Short Summary
    -------------
    Calculate the number of rows per data section to process.

    Parameters
    ----------
    model: instance of Data Model
       DM object for input

    buffsize: int
       size of data section (buffer) in bytes

    cubeshape: (int, int, int) tuple
       shape of input dataset

    nreads: int
       number of reads in input dataset

    Returns
    -------
    nrows: int
       number of rows in buffer of data section
    """
    bitpix = model.data.dtype.itemsize
    bytepix = int(abs(bitpix) / 8)
    if bytepix < 1:
        bytepix = 1

    nrows = int(buffsize / (bytepix * cubeshape[2] * nreads))
    if nrows < 1:
        nrows = 1
    if nrows > cubeshape[1]:
        nrows = cubeshape[1]

    return nrows


def calc_slope(data_sect, gdq_sect, frame_time, opt_res, save_opt, rn_sect,
               gain_sect, i_max_seg, ngroups, weighting, f_max_seg, group_time):
    """
    Short Summary
    -------------
    Compute the slope of each segment for each pixel in the data cube section
    for the current integration. Each segment has its slope fit in fit_lines();
    that slope and other quantities from the fit are added to the 'optional
    result' object by append_arr() from the appropriate 'CASE' (type of segment)
    in fit_next_segment().

    Parameters
    ----------
    data_sect: float
        section of input data cube array

    gdq_sect: float
        section of GROUPDQ data quality array

    frame_time: float
        integration time

    opt_res: OptRes object
        contains quantities related to fitting for optional output

    save_opt: boolean
       save optional fitting results

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    i_max_seg: int
        used for size of initial allocation of arrays for optional results;
        maximum possible number of segments within the ramp, based on the
        number of CR flags

    ngroups: int
        number of groups per integration

    weighting: string
        'optimal' specifies that optimal weighting should be used; currently
        the only weighting supported.

    f_max_seg: int
        actual maximum number of segments within a ramp, based on the fitting
        of all ramps; later used when truncating arrays before output.

    group_time: float
        Time increment between groups, in seconds.

    Returns
    -------
    gdq_sect: int, 3D array
        data quality flags for pixels in section

    inv_var: float, 1D array
        values of 1/variance for good pixels

    opt_res: OptRes object
        contains quantities related to fitting for optional output

    f_max_seg: int
        actual maximum number of segments within a ramp, updated here based on
        fitting ramps in the current data section; later used when truncating
        arrays before output.

    num_seg: int, 1D array
        numbers of segments for good pixels
    """
    nreads, asize2, asize1 = data_sect.shape
    npix = asize2 * asize1  # number of pixels in section of 2D array

    all_pix = np.arange(npix)
    arange_nreads_col = np.arange(nreads)[:, np.newaxis]
    start = np.zeros(npix, dtype=np.int32) # lowest channel in fit

    # Highest channel in fit initialized to last read
    end = np.zeros(npix, dtype=np.int32) + (nreads - 1)

    pixel_done = (end < 0) # False until processing is done

    inv_var = np.zeros(npix, dtype=np.float32) # inverse of fit variance
    num_seg = np.zeros(npix, dtype=np.int32) # number of segments per pixel

    # End stack array - endpoints for each pixel
    # initialize with nreads for each pixel; set 1st channel to 0
    end_st = np.zeros((nreads + 1, npix), dtype=np.int32)
    end_st[0, :] = nreads - 1

    # end_heads is initially a tuple populated with every pixel that is
    # either saturated or contains a cosmic ray based on the input DQ
    # array, so is sized to accomodate the maximum possible number of
    # pixels flagged. It is later compressed to be an array denoting
    # the number of endpoints per pixel.
    end_heads = np.ones(npix * nreads, dtype=np.int32)

    # Create nominal 2D ERR array, which is 1st slice of
    #    avged_data_cube * readtime
    err_2d_array = data_sect[0, :, :] * frame_time
    err_2d_array[err_2d_array < 0] = 0

    # Frames >= start and <= end will be masked. However, the first channel
    #   to be included in fit will be the read in which a cosmic ray has
    #   been flagged
    mask_2d = ((arange_nreads_col >= start[np.newaxis, :]) &
               (arange_nreads_col <= end[np.newaxis, :]))

    end = 0 # array no longer needed

    # Section of GROUPDQ dq section, excluding bad dq values in mask
    gdq_sect_r = np.reshape(gdq_sect, (nreads, npix))
    mask_2d[gdq_sect_r != 0] = False  # saturated or CR-affected
    mask_2d_init = mask_2d.copy() # initial flags for entire ramp

    wh_f = np.where(np.logical_not(mask_2d))

    these_p = wh_f[1] # coordinates of pixels flagged as False
    these_r = wh_f[0] # reads of pixels flagged as False

    del wh_f

    # Populate end_st to contain the set of end points for each pixel.
    # Populate end_heads to initially include every pixel that is either
    # saturated or contains a cosmic ray. Skips the duplicated final group
    # for saturated pixels. Saturated pixels resulting in a contiguous set
    # of intervals of length 1 will later be flagged as too short
    # to fit well.
    for ii, val in enumerate(these_p):
        if (these_r[ii] != (nreads - 1)):
            end_st[end_heads[these_p[ii]], these_p[ii]] = these_r[ii]
            end_heads[these_p[ii]] += 1

    # Sort and reverse array to handle the order that saturated pixels
    # were added
    end_st.sort(axis=0)
    end_st = end_st[::-1]

    # Reformat to designate the number of endpoints per pixel; compress
    # to specify number of groups per pixel
    end_heads = (end_st > 0).sum(axis=0)

    # Create object to hold optional results
    opt_res.init_2d(npix, i_max_seg, save_opt)

    # LS fit until 'nreads' iterations or all pixels in
    #    section have been processed
    for iter_num in range(nreads):
        if pixel_done.all():
            break

       # frames >= start and <= end_st will be included in fit
        mask_2d = ((arange_nreads_col >= start) &
                   (arange_nreads_col <
                    (end_st[end_heads[all_pix] - 1, all_pix] + 1)))

        mask_2d[gdq_sect_r != 0] = False # RE-exclude bad group dq values

        # for all pixels, update arrays, summing slope and variance
        f_max_seg, num_seg = \
              fit_next_segment(start, end_st, end_heads, pixel_done, data_sect,
                 mask_2d, mask_2d_init, inv_var, num_seg, opt_res, save_opt,
                 rn_sect, gain_sect, ngroups, weighting, f_max_seg)

        if f_max_seg is None:
            f_max_seg = 1

    arange_nreads_col = 0
    all_pix = 0

    return gdq_sect, inv_var, opt_res, f_max_seg, num_seg


def fit_next_segment(start, end_st, end_heads, pixel_done, data_sect, mask_2d,
                     mask_2d_init, inv_var, num_seg, opt_res, save_opt, rn_sect, 
                     gain_sect, ngroups, weighting, f_max_seg):
    """
    Extended Summary
    ----------------
    Call routine to LS fit masked data for a single segment for all pixels in
    data section. Then categorize each pixel's fitting interval based on
    interval length, and whether the interval is at the end of the array.
    Update the start array, the end stack array, the end_heads array which
    contains the number of endpoints. For pixels in which the fitting intervals
    are long enough, the resulting slope and variance are added to the
    appropriate stack arrays.  The first channel to fit in a segment is either
    the first group in the ramp, or a group in which a cosmic ray has been
    flagged.

    Parameters
    ----------
    start: int, 1D array
        lowest channel in fit

    end_st: int, 2D array
        stack array of endpoints

    end_heads: int, 1D array
        number of endpoints for each pixel

    pixel_done: boolean, 1D array
        whether each pixel's calculations are completed

    data_sect: float, 3D array
        data cube section

    mask_2d: bool, 2D array
        delineates which channels to fit for each pixel

    mask_2d_init: bool, 2D array
        copy of intial mask_2d

    inv_var: float, 1D array
        values of 1/variance for good pixels

    num_seg: int, 1D array
        numbers of segments for good pixels

    opt_res: OptRes object
        optional fitting results to output

    save_opt: boolean
       save optional fitting results

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    ngroups: int
        number of groups per integration

    weighting: string
        'optimal' specifies that optimal weighting should be used; currently
        the only weighting supported.

    f_max_seg: int
        actual maximum number of segments within a ramp, updated here based on
        fitting ramps in the current data section; later used when truncating
        arrays before output.

    Returns
    -------
    f_max_seg: int
        actual maximum number of segments within a ramp, updated here based on
        fitting ramps in the current data section; later used when truncating
        arrays before output.

    num_seg: int, 1D array
        numbers of segments for good pixels
    """
    nreads, asize2, asize1 = data_sect.shape # Note: nreads is a scalar here
    all_pix = np.arange(asize2 * asize1)

    ramp_mask_sum = mask_2d_init.sum(axis=0)

    # Each returned array below is 1D, for all npix pixels for current segment
    slope, intercept, variance, sig_intercept, sig_slope = fit_lines(data_sect,
              mask_2d, rn_sect, gain_sect, ngroups, weighting)

    end_locs = end_st[end_heads[all_pix] - 1, all_pix]
    l_interval = end_locs - start # fitting interval length

    wh_done = (start == -1) # done pixels
    l_interval[wh_done] = 0  # set interval lengths for done pixels to 0

    # Create array to set when each good pixel is classified for the current
    #   semiramp (to enable unclassified pixels to have their arrays updated)
    got_case = np.zeros((asize1*asize2), dtype=np.bool)

    # CASE A) Long enough (semiramp has >2 groups), at end of ramp
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - add slopes and variances to running sums
    wh_check = np.where((l_interval>2) & (end_locs==nreads - 1) & (~pixel_done))

    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        start[these_pix] = -1   # all processing for this pixel is completed
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] = 0
        pixel_done[these_pix] = True # all processing for pixel is completed
        got_case[ these_pix ] = True

        g_pix = these_pix[variance[these_pix] > 0.] # good pixels
        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]

            # Append results to arrays
            opt_res.append_arr(num_seg, g_pix, intercept, slope,
                sig_intercept, sig_slope, inv_var, save_opt)

            num_seg[g_pix] += 1
            f_max_seg = max(f_max_seg, num_seg.max())


    # CASE B) Long enough (semiramp has >2 ), not at array end (meaning final
    #    group for this semiramp is not final group of the whole ramp)
    #    - remove current end from end stack
    #    - decrement number of ends
    #    - add slopes and variances to running sums
    wh_check = np.where((l_interval > 2) & (end_locs != nreads - 1) & ~pixel_done)
    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        got_case[ these_pix ] = True

        start[these_pix] = end_locs[these_pix]
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] -= 1
        end_heads[end_heads < 0.] = 0.

        g_pix = these_pix[variance[these_pix] > 0.] # good pixels

        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]

            # Append results to arrays
            opt_res.append_arr(num_seg, g_pix, intercept, slope,
                sig_intercept, sig_slope, inv_var, save_opt)

            num_seg[g_pix] += 1
            f_max_seg = max(f_max_seg, num_seg.max())

            # If there are pixels with no later good groups, update stack
            #   arrays accordingly
            c_mask_2d_init = mask_2d_init.copy()

            # create array: 0...nreads-1 in a column for each pixel
            arr_ind_all = np.array( [np.arange(nreads),] *
                                    c_mask_2d_init.shape[1]).transpose()
            wh_c_start_all = np.zeros( c_mask_2d_init.shape[1], dtype=np.uint8)
            wh_c_start_all[ g_pix ]= start[ g_pix]

            # set to False all groups before start group
            c_mask_2d_init[ arr_ind_all < wh_c_start_all ] = False

            # select pixels having all groups False from start to ramp end
            wh_rest_false = np.where( c_mask_2d_init.sum(axis=0) == 0)
            if(len(wh_rest_false[0]) > 0):
                pix_rest_false = wh_rest_false[0]
                start[ pix_rest_false ] = -1
                end_st[ end_heads[ pix_rest_false ] - 1, pix_rest_false ] = 0
                end_heads[ pix_rest_false ] = 0
                pixel_done[ pix_rest_false ] = True # all processing is complete


    # CASE C) - dataset has NGROUPS=1 ; so special fitting is done for all pixels
    #    and all intervals are at the end of the array.
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - add slopes and variances to running sums
    #    - set pixel_done to True to designate all fitting done
    if (ngroups == 1):
        start[all_pix] = -1
        end_st[end_heads[all_pix] - 1, all_pix] = 0
        end_heads[all_pix] = 0
        pixel_done[all_pix] = True

        wh_check = np.where(mask_2d_init[0, :] & (ramp_mask_sum == 1))
        if(len(wh_check[0]) > 0):
            g_pix = wh_check[0]

        # Ignore all pixels having no good groups (so the single group is bad)
        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]

            # Append results to arrays
            opt_res.append_arr(num_seg, g_pix, intercept, slope,
                sig_intercept, sig_slope, inv_var, save_opt)

            num_seg[g_pix] = 1

        return 1, num_seg


    # CASE D) - dataset has NGROUPS=2, so special fitting is done for all pixels.
    #    All segments are at the end of the array.
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - add slopes and variances to running sums
    #    - set pixel_done to True to designate all fitting done
    if (ngroups == 2):
        start[all_pix] = -1
        end_st[end_heads[all_pix] - 1, all_pix] = 0
        end_heads[all_pix] = 0
        pixel_done[all_pix] = True

        g_pix = all_pix[variance[all_pix] > 0.]
        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]

            opt_res.append_arr(num_seg, g_pix, intercept, slope,
                sig_intercept, sig_slope, inv_var, save_opt)

            num_seg[g_pix] = 1

            return 1, num_seg


    # CASE E) - interval too short to fit normally (only 2 good groups)
    #    At end of array, NGROUPS>1, but exclude NGROUPS==2 datasets
    #    as they are covered in CASE D
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - add slopes and variances to running sums
    #    - set pixel_done to True to designate all fitting done
    wh_check = np.where((l_interval == 1) & (end_locs == nreads - 1) &
                        (nreads > 1) & (ngroups != 2) & (~pixel_done))

    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        got_case[ these_pix ] = True

        start[these_pix] = -1
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] = 0
        pixel_done[these_pix] = True

        g_pix = these_pix[variance[these_pix] > 0.] # good pixels

        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]

            # Append results to arrays
            opt_res.append_arr(num_seg, g_pix, intercept, slope,
                               sig_intercept, sig_slope, inv_var, save_opt)

            num_seg[g_pix] += 1
            f_max_seg = max(f_max_seg, num_seg.max())


    # CASE F) - full-length ramp has 2 good groups not at array end
    #    - use the 2 good reads to get the slope
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of end to 0
    #    - add slopes and variances to running sums
    #    - set pixel_done to True to designate all fitting done

    # Copy mask, as will modify when calculating the number of later good groups
    c_mask_2d_init = mask_2d_init.copy()

    wh_check = np.where((l_interval == 2) & ( ngroups >2 ) &
                        (end_locs != nreads - 1) & ~pixel_done)
    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        got_case[ these_pix ] = True
        inv_var[these_pix] += 1.0 / variance[these_pix]

        # create array: 0...nreads-1 in a column for each pixel
        arr_ind_all = np.array([np.arange(nreads),] *
                               c_mask_2d_init.shape[1]).transpose()
        wh_c_start_all = np.zeros( mask_2d_init.shape[1], dtype=np.uint8)
        wh_c_start_all[ these_pix ]= start[ these_pix]

        # set to False all groups before start group
        c_mask_2d_init[ arr_ind_all < wh_c_start_all ] = 0

        tot_good_groups = c_mask_2d_init.sum(axis=0 )

        # Select pixels having at least 2 later good groups (these later good
        #   groups are a segment whose slope will be calculated)
        wh_more = np.where( tot_good_groups[these_pix] > 1 )
        pix_more = these_pix[ wh_more ]
        start[ pix_more ] = end_locs[ pix_more ]
        end_st[ end_heads[ pix_more ] - 1, pix_more ] = 0
        end_heads[ pix_more ] -= 1

        # Select pixels having less than 2 later good groups (these later good
        #   groups will not be used)
        wh_only = np.where( tot_good_groups[these_pix] <= 1 )
        pix_only = these_pix[ wh_only ]
        start[ pix_only ] = -1
        end_st[ end_heads[ pix_only ] - 1, pix_only ] = 0
        end_heads[ pix_only ] = 0
        pixel_done[ pix_only ] = True # all processing for pixel is completed

        end_heads[(end_heads < 0.)] = 0.

        # Append results to arrays
        opt_res.append_arr(num_seg, these_pix, intercept, slope,
                           sig_intercept, sig_slope, inv_var, save_opt)

        num_seg[these_pix] += 1
        f_max_seg = max(f_max_seg, num_seg.max())


    # CASE G) - full-length ramp has a good group on 0th group of the entire ramp,
    #    and no later good groups. Will use single good group data as the slope
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of end to 0
    #    - add slopes and variances to running sums
    #    - set pixel_done to True to designate all fitting done
    wh_check = np.where(mask_2d_init[0, :] & ~mask_2d_init[1, :] &
                   (ramp_mask_sum == 1) & ~pixel_done)

    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        got_case[ these_pix ] = True

        start[these_pix] = -1
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] = 0
        pixel_done[these_pix] = True # all processing for pixel is completed
        inv_var[these_pix] += 1.0 / variance[these_pix]

        # Append results to arrays
        opt_res.append_arr(num_seg, these_pix, intercept, slope,
                           sig_intercept, sig_slope, inv_var, save_opt)

        num_seg[these_pix] += 1
        f_max_seg = max(f_max_seg, num_seg.max())

    # CASE H) - the segment has a good 0th group and a bad 1st group. For the
    #  data from the 0th good group of this segment to possibly be used as a
    #  slope, that group must necessarily be the 0th group of the entire ramp.
    #  It is possible to have a single 'good' group segment after the 0th group
    #  of the ramp; in that case the 0th group and the 1st group would both have
    #  to be CRs, and the data of the 0th group would not be included as a slope.
    #  For a good 0th group in a ramp followed by a bad 1st group there must be
    #  good groups later in the segment because if there were not, the segment
    #  would already have be classified as CASE G. In this situation, since
    #  there are later good groups in the segment, those later good groups will
    #  be used in the slope computation, and the 0th good group will not be.
    #  As a result, for all instances of these types of segments, the data in the
    #  initial good group will not be used in the slope calculation, but the
    #  arrays for the indices for the ramp (end_st, etc) are appropriately
    #  adjusted.
    #   - increment start array
    #   - remove current end from end stack
    #   - decrement number of ends
    wh_check = np.where(mask_2d_init[0, :] & ~mask_2d_init[1, :] & ~pixel_done)

    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        got_case[ these_pix ] = True
        start[ these_pix ] += 1
        start[ start > nreads-1 ] = nreads -1  # to keep at max level
        end_st[ end_heads[ these_pix ] - 1, these_pix ] = 0
        end_heads[ these_pix ] -= 1
        end_heads [end_heads < 0. ] = 0.


    # CASE OTHER) - all other types of segments not covered earlier. No segments
    #   handled here have adequate data, but the stack arrays are updated.
    #   - increment start array
    #   - remove current end from end stack
    #   - decrement number of ends
    wh_check = np.asarray( np.where( ~pixel_done & ~got_case ))

    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        start[ these_pix ] += 1
        start[ start > nreads-1 ] = nreads -1  # to keep at max level
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] -= 1
        end_heads[end_heads < 0.] = 0.

    return f_max_seg, num_seg


def fit_lines(data, mask_2d, rn_sect, gain_sect, ngroups, weighting):
    """
    Extended Summary
    ----------------
    Do linear least squares fit to data cube in this integration for a single
    segment for all pixels.  In addition to applying the mask due to identified
    cosmic rays, the data is also masked to exclude intervals that are too short
    to fit well. The first channel to fit in a segment is either the first group
    in the ramp, or a group in which a cosmic ray has been flagged.

    Parameters
    ----------
    data: float
       array of values for current data section

    mask_2d: boolean, 2D array
       delineates which channels to fit for each pixel

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    ngroups: int
        number of groups per integration

    weighting: string
        'optimal' specifies that optimal weighting should be used; currently
        the only weighting supported.

    Returns
    -------
    Note: all of these pertain to a single segment (hence '_s')

    slope_s: float, 1D array
       weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
       y-intercepts from fit for data section

    variance_s: float, 1D array
       variance of residuals for fit for data section

    sig_intercept_s: float, 1D array
       sigma of y-intercepts from fit for data section

    sig_slope_s: float, 1D array
       sigma of slopes from fit for data section (for a single segment)

    """
    # verify that incoming data is either 2 or 3-dimensional
    try:
        assert (data.ndim == 2 or data.ndim == 3)
    except AssertionError:
        log.error('FATAL ERROR: Data input to fit_lines must be 2 or 3 dimensions')

    # To ensure that the first channel to be fit is the cosmic-ray-affected
    #   group, the channel previous to each channel masked as good is
    #   also masked as good. This is only for the local purpose of setting
    #   the first channel, and will not propagate beyond this current function
    #   call.
    c_mask_2d = mask_2d.copy()
    wh_mask_2d = np.where(c_mask_2d)
    c_mask_2d[np.maximum(wh_mask_2d[0] - 1, 0), wh_mask_2d[1]] = True

    del wh_mask_2d

    # num of reads/pixel unmasked
    nreads_1d = c_mask_2d.astype(np.int16).sum(axis=0)
    npix = c_mask_2d.shape[1]

    slope_s = np.zeros(npix, dtype=np.float32)
    variance_s = np.zeros(npix, dtype=np.float32)
    intercept_s = np.zeros(npix, dtype=np.float32)
    sig_intercept_s = np.zeros(npix, dtype=np.float32)
    sig_slope_s = np.zeros(npix, dtype=np.float32)

    # Calculate slopes etc. for datasets having either 1 or 2 groups per
    #   integration, and return
    if (ngroups == 1): # process all pixels in 1 group/integration dataset
        slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s = \
           fit_1_group(slope_s, intercept_s, variance_s, sig_intercept_s,
                       sig_slope_s, npix, data, c_mask_2d)

        return slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s

    if (ngroups == 2): # process all pixels in 2 group/integration dataset
        rn_sect_1d = rn_sect.reshape(npix)
        slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s = \
           fit_2_group(slope_s, intercept_s, variance_s, sig_intercept_s,
                       sig_slope_s, npix, data, c_mask_2d, rn_sect_1d)

        return slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s

    # reshape data_masked
    data_masked = data * np.reshape(c_mask_2d, data.shape)
    data_masked = np.reshape(data_masked, (data_masked.shape[0], npix))

   # For datasets having >2 groups/integration, for any semiramp in which the
   #   0th group is good and the 1st group is bad, determine whether or not to
   #   use the 0th group.
    wh_pix_1r = np.where(c_mask_2d[0,:] & (np.logical_not(c_mask_2d[1,:])))

    if (len(wh_pix_1r[0]) > 0 ):
        slope_s, intercept_s, variance_s, sig_intercept_s, \
            sig_slope_s = fit_single_read(slope_s, intercept_s,
            variance_s, sig_intercept_s, sig_slope_s, npix, data, wh_pix_1r)

    del wh_pix_1r

   # For datasets having >2 groups/integrations, for any semiramp in which only
   #   the 0th and 1st group are good, set slope, etc
    wh_pix_2r = np.where( c_mask_2d.sum(axis=0) ==2) # ramps with 2 good groups
    slope_s, intercept_s, variance_s, sig_slope_s, sig_intercept_s = \
        fit_double_read( c_mask_2d, wh_pix_2r, data_masked, slope_s, intercept_s,
        variance_s, sig_slope_s, sig_intercept_s, rn_sect)

    del wh_pix_2r

    # Select ramps having >2 good groups
    wh_pix_to_use = np.where(c_mask_2d.sum(axis=0) > 2)

    good_pix = wh_pix_to_use[0]  # Ramps with >2 good groups
    data_masked = data_masked[:, good_pix]

    del wh_pix_to_use

    xvalues = np.arange(data_masked.shape[0])[:, np.newaxis] * c_mask_2d
    xvalues = xvalues[:, good_pix]  # set to those pixels to be used

    c_mask_2d = c_mask_2d[:, good_pix]
    nreads_1d = nreads_1d[good_pix]

    if weighting.lower() == 'optimal': # fit using optimal weighting
        # get sums from optimal weighting
        sumx, sumxx, sumxy, sumy, nreads_wtd, xvalues = calc_opt_sums( rn_sect,
                gain_sect, data_masked, c_mask_2d, xvalues, good_pix )

        slope, intercept, sig_slope, sig_intercept = calc_opt_fit( nreads_wtd,
                      sumxx, sumx, sumxy, sumy)

        variance = sig_slope**2.  # variance due to fit values

    elif weighting.lower() == 'unweighted': # fit using unweighted weighting
        # get sums from unweighted weighting
        sumx, sumxx, sumxy, sumy =\
              calc_unwtd_sums(data_masked, xvalues)

        slope, intercept, sig_slope, sig_intercept, line_fit =\
               calc_unwtd_fit(xvalues, nreads_1d, sumxx, sumx, sumxy, sumy)

        denominator = nreads_1d * sumxx - sumx**2

        variance = nreads_1d / denominator
        denominator = 0

    else: # unsupported weighting type specified
        log.error('FATAL ERROR: unsupported weighting type specified.')

    line_fit = 0

    slope_s[good_pix] = slope
    variance_s[good_pix] = variance
    intercept_s[good_pix] = intercept
    sig_intercept_s[good_pix] = sig_intercept
    sig_slope_s[good_pix] = sig_slope

    return slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s


def fit_single_read(slope_s, intercept_s, variance_s, sig_intercept_s,
                    sig_slope_s, npix, data, wh_pix_1r):
    """
    Short Summary
    -------------
    For datasets having >2 groups/integrations, for any semiramp in which the
    0th group is good and the 1st group is either SAT or CR, set slope, etc.

    Parameters
    ----------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section

    npix: int
        number of pixels in 2D array

    data: float
        array of values for current data section

    wh_pix_1r: tuple
        locations of pixels whose only good group is the 0th group

    Returns
    -------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section
    """
    data0_slice = data[0, :, :].reshape(npix)
    slope_s[wh_pix_1r] = data0_slice[wh_pix_1r]

    # The following arrays will have values correctly calculated later; for
    #    now they are just place-holders
    variance_s[wh_pix_1r] = utils.LARGE_VARIANCE
    sig_slope_s[wh_pix_1r] = 0.
    intercept_s[wh_pix_1r] = 0.
    sig_intercept_s[wh_pix_1r] = 0.

    return slope_s, intercept_s, variance_s, sig_slope_s, sig_intercept_s


def fit_double_read(mask_2d, wh_pix_2r, data_masked, slope_s, intercept_s,
                    variance_s, sig_slope_s, sig_intercept_s, rn_sect):
    """
    Short Summary
    -------------
    Process all semi-ramps having exactly 2 good groups. May need to optimize
    later to remove loop over pixels.

    Parameters
    ----------
    mask_2d: bool, 2D array
        delineates which channels to fit for each pixel

    wh_pix_2r: tuple
        locations of pixels whose only good groups are the 0th and the 1st

    data_masked: float, 2D array
        masked values for all pixels in data section

    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    Returns
    -------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section
    """
    for ff in range(len(wh_pix_2r[0])): # loop over the pixels
        pixel_ff = wh_pix_2r[0][ff] # pixel index (1d)
        rn = rn_sect.flatten()[pixel_ff] # read noise for this pixel

        read_nums = np.where( mask_2d[:,pixel_ff])
        second_read = read_nums[0][1]
        data_ramp = data_masked[:, pixel_ff] * mask_2d[:, pixel_ff]
        data_semi = data_ramp[ mask_2d[:, pixel_ff]] # picks only the 2
        diff_data =  data_semi[1] - data_semi[0]

        slope_s[ pixel_ff ] = diff_data
        intercept_s[ pixel_ff ] = data_semi[1]*(1.- second_read) + \
            data_semi[0]*second_read # by geometry
        variance_s[ pixel_ff ] = 2.0 * rn * rn
        sig_slope_s[pixel_ff] = np.sqrt(2) * rn
        sig_intercept_s[ pixel_ff ] = np.sqrt(2) * rn

    return slope_s, intercept_s, variance_s, sig_slope_s, sig_intercept_s


def calc_unwtd_fit(xvalues, nreads_1d, sumxx, sumx, sumxy, sumy):
    """
    Extended Summary
    ----------------
    Do linear least squares fit to data cube in this integration, using
    unweighted fits to the segments. Currently not supported.

    Parameters
    ----------
    xvalues: int, 1D array
        indices of valid pixel values for all groups

    nreads_1d: int, 1D array
        number of reads in an integration

    sumxx: float
        sum of squares of xvalues

    sumx: float
        sum of xvalues

    sumxy: float
        sum of product of xvalues and data

    sumy: float
        sum of data

    Returns
    -------
    slope: float, 1D array
       weighted slope for current iteration's pixels for data section

    intercept: float, 1D array
       y-intercepts from fit for data section

    sig_slope: float, 1D array
       sigma of slopes from fit for data section

    sig_intercept: float, 1D array
       sigma of y-intercepts from fit for data section

    line_fit: float, 1D array
       values of fit using slope and intercept
    """

    denominator = nreads_1d * sumxx - sumx**2
    slope = (nreads_1d * sumxy - sumx * sumy) / denominator
    intercept = (sumxx * sumy - sumx * sumxy) / denominator
    sig_intercept = (sumxx / denominator)**0.5
    sig_slope = (nreads_1d / denominator)**0.5

    line_fit = (slope * xvalues) + intercept

    return slope, intercept, sig_slope, sig_intercept, line_fit


def calc_opt_fit(nreads_wtd, sumxx, sumx, sumxy, sumy):
    """
    Extended Summary
    ----------------
    Do linear least squares fit to data cube in this integration for a single
    semi-ramp for all pixels, using optimally weighted fits to the semi_ramps.
    The weighting uses the formulation by Fixsen (Fixsen et al, PASP, 112, 1350).
    Note - these weights, sigmas, and variances pertain only to the fitting, and
    the variances are *NOT* the variances of the slope due to noise.

    Parameters
    ----------
    nreads_wtd: float, 1D array
        sum of product of data and optimal weight

    sumxx: float
        sum of squares of xvalues

    sumx: float
        sum of xvalues

    sumxy: float
        sum of product of xvalues and data

    sumy: float
        sum of data

    Returns
    -------
    slope: float, 1D array
       weighted slope for current iteration's pixels for data section

    intercept: float, 1D array
       y-intercepts from fit for data section

    sig_slope: float, 1D array
       sigma of slopes from fit for data section

    sig_intercept: float, 1D array
       sigma of y-intercepts from fit for data section
    """
    denominator = nreads_wtd * sumxx - sumx**2
    slope = (nreads_wtd * sumxy - sumx * sumy) / denominator
    intercept = (sumxx * sumy - sumx * sumxy) / denominator
    sig_intercept = (sumxx / denominator)**0.5
    sig_slope = (nreads_wtd / denominator)**0.5 # STD of the slope's fit

    return slope, intercept, sig_slope, sig_intercept


def fit_1_group(slope_s, intercept_s, variance_s, sig_intercept_s,
                sig_slope_s, npix, data, mask_2d):
    """
    Extended Summary
    ----------------
    This function sets the fitting arrays for datasets having only 1 group
    per integration.

    Parameters
    ----------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section

    npix: int
        number of pixels in 2d array

    data: float
        array of values for current data section

    mask_2d: bool, 2D array
        delineates which channels to fit for each pixel

    Returns
    -------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section
    """
    # For pixels not saturated, recalculate the slope as the value of the SCI
    # data in that group, which will later be divided by the group exposure
    # time to give the count rate. Recalculate other fit quantities to be
    # benign.
    slope_s = data[0, :, :].reshape(npix)

    # The following arrays will have values correctly calculated later; for
    #    now they are just place-holders
    variance_s = np.zeros(npix, dtype=np.float32) + utils.LARGE_VARIANCE
    sig_slope_s = slope_s * 0.
    intercept_s = slope_s * 0.
    sig_intercept_s = slope_s * 0.

    # For saturated pixels, overwrite slope with benign values.
    wh_sat0 = np.where(np.logical_not(mask_2d[0, :]))

    if (len(wh_sat0[0]) > 0):
        sat_pix = wh_sat0[0]
        slope_s[sat_pix] = 0.

    return slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s


def fit_2_group(slope_s, intercept_s, variance_s, sig_intercept_s,
                sig_slope_s, npix, data, mask_2d, rn_sect_1d):
    """
    Extended Summary
    ----------------
    This function sets the fitting arrays for datasets having only 2 groups
    per integration.

    Parameters
    ----------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section

    npix: int
        number of pixels in 2d array

    data: float
        array of values for current data section

    mask_2d: bool, 2D array
        delineates which channels to fit for each pixel

    rn_sect_1d: float, 1D array
        read noise values for all pixels in data section

    Returns
    -------
    slope_s: float, 1D array
        weighted slope for current iteration's pixels for data section

    intercept_s: float, 1D array
        y-intercepts from fit for data section

    variance_s: float, 1D array
        variance of residuals for fit for data section

    sig_intercept_s: float, 1D array
        sigma of y-intercepts from fit for data section

    sig_slope_s: float, 1D array
        sigma of slopes from fit for data section
    """
    # For pixels saturated on the first group, overwrite fit values with
    # benign values to be recalculated later.
    wh_sat0 = np.where(np.logical_not(mask_2d[0, :]))

    if (len(wh_sat0[0]) > 0):
        sat_pix = wh_sat0[0]
        slope_s[sat_pix] = 0.
        variance_s[sat_pix] = 0.
        sig_slope_s[sat_pix] = 0.
        intercept_s[sat_pix] = 0.
        sig_intercept_s[sat_pix] = 0.
    del wh_sat0

    # For pixels saturated on the second group, recalculate the slope as
    # the value of the SCI data in the first group, which will later be
    # divided by the group exposure time to give the count rate, and
    # recalculate the other fit quantities to be benign. Note: these pixels
    # will already have been handled earlier (for intervals of arbitrary
    # length) in this function, but are being included here to explicitly
    # cover all possibilities for pixels in datasets with ngroups=2. Will
    # later consider refactoring.
    wh_sat1 = np.where((mask_2d[:, :].sum(axis=0) == 1) & mask_2d[0, :])

    if (len(wh_sat1[0]) > 0):
        data0_slice = data[0, :, :].reshape(npix)
        slope_s[wh_sat1] = data0_slice[wh_sat1]
        variance_s[wh_sat1] = 0.
        sig_slope_s[wh_sat1] = 0.
        intercept_s[wh_sat1] = 0.
        sig_intercept_s[wh_sat1] = 0.
    del wh_sat1

    # For pixels with no saturated values, recalculate the slope as the
    # difference between the values of the second and first groups (1-based),
    # which will later be divided by the group exposure time to give the count
    # rate, and recalculate other fit quantities to be benign.
    wh_sat_no = np.where(mask_2d[:, :].sum(axis=0) == 2)

    if (len(wh_sat_no[0]) > 0):
        data0_slice = data[0, :, :].reshape(npix)
        data1_slice = data[1, :, :].reshape(npix)
        slope_s[wh_sat_no] = data1_slice[wh_sat_no] - data0_slice[wh_sat_no]
        sig_slope_s[wh_sat_no] = np.sqrt(2) * rn_sect_1d[wh_sat_no]
        intercept_s[wh_sat_no] = data0_slice[wh_sat_no] -\
                data1_slice[wh_sat_no] # by geometry
        sig_intercept_s[wh_sat_no] = np.sqrt(2) * rn_sect_1d[wh_sat_no]
        variance_s[wh_sat_no] = np.sqrt(2) * rn_sect_1d[wh_sat_no]

    del wh_sat_no

    return slope_s, intercept_s, variance_s, sig_intercept_s, sig_slope_s


def calc_num_seg(gdq, n_int):
    """
    Extended Summary
    ----------------
    Calculate the maximum number of segments that will be be fit within an
    integration, calculated over all pixels and all integrations.  This value
    is based on the locations of cosmic ray-affected pixels in all of the ramps,
    and will be used to allocate arrays used for the optional output product.

    Parameters
    ----------
    gdq: float, 3D array
        cube of GROUPDQ array for a data

    n_int: int
        total number of integrations in data set

    Return:
    -------
    int(max_cr) +1; int
        maxmimum number of segments; n CRS implies n+1 segments
    """
    max_cr = 0  # max number of CRS for all integrations

    # For all 2d pixels, get max number of CRs along their ramps
    # to use as a surrogate for the number of semiramps along the ramps
    for nint in range(n_int):  # loop over integrations
        gdq_cr = np.bitwise_and( gdq[nint, :, :, :], dqflags.group['JUMP_DET'])
        temp_max_cr = int((gdq_cr.sum(axis=0)).max()/dqflags.group['JUMP_DET'])
        max_cr = max( max_cr, temp_max_cr )

    return int(max_cr) +1 # n CRS implies n+1 segments


def calc_unwtd_sums(data_masked, xvalues):
    """
    Short Summary
    -------------
    Calculate the sums needed to determine the slope and intercept (and sigma
    of each) using an unweighted fit. Unweighted fitting currently not
    supported.

    Parameters
    ----------
    data_masked: float, 2D array
        masked values for all pixels in data section

    xvalues: int, 1D array
        indices of valid pixel values for all groups

    Return:
    -------
    sumx: float
        sum of xvalues

    sumxx: float
        sum of squares of xvalues

    sumxy: float
        sum of product of xvalues and data

    sumy: float
        sum of data

    """
    sumx = xvalues.sum(axis=0)
    sumxx = (xvalues**2).sum(axis=0)
    sumy = (np.reshape(data_masked.sum(axis=0), sumx.shape))
    sumxy = (xvalues * np.reshape(data_masked, xvalues.shape)).sum(axis=0)

    return sumx, sumxx, sumxy, sumy


def calc_opt_sums(rn_sect, gain_sect, data_masked, mask_2d, xvalues, good_pix):
    """
    Short Summary
    -------------
    Calculate the sums needed to determine the slope and intercept (and sigma of
    each) using the optimal weights.  For each good pixel's segment, from the
    initial and final indices and the corresponding number of counts, calculate
    the SNR. From the SNR, calculate the weighting exponent using the formulation
    by Fixsen (Fixsen et al, PASP, 112, 1350). Using this exponent and the gain
    and the readnoise, the weights are calculated from which the sums are
    calculated.

    Parameters
    ----------
    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    data_masked: float, 2D array
        masked values for all pixels in data section

    mask_2d: bool, 2D array
        delineates which channels to fit for each pixel

    xvalues: int, 2D array
        indices of valid pixel values for all groups

    good_pix: int, 1D array
        indices of pixels having valid data for all groups

    Return:
    -------
    sumx: float
        sum of xvalues

    sumxx: float
        sum of squares of xvalues

    sumxy: float
        sum of product of xvalues and data

    sumy: float
        sum of data

    nreads_wtd: float, 1D array
        sum of optimal weights

    xvalues: int, 2D array
        rolled up indices of valid pixel values for all groups
    """
    c_mask_2d = mask_2d.copy() # copy the mask to prevent propagation

    # Return 'empty' sums if there is no more data to fit
    if (data_masked.size == 0):
        return np.array([]), np.array([]), np.array([]), np.array([]),\
               np.array([]), np.array([])

    # get initial group for each good pixel for this semiramp
    fnz = np.argmax(c_mask_2d, axis=0)

    # For those pixels that are all False, set to sentinel value of -1
    fnz[c_mask_2d.sum(axis=0) == 0] = -1

    mask_2d_sum = c_mask_2d.sum(axis=0)   # number of valid groups/pixel

    # get final valid group for each pixel for this semiramp
    ind_lastnz = fnz + mask_2d_sum - 1

    # get SCI value of initial good group for semiramp
    data_zero = data_masked[fnz, range(data_masked.shape[1])]

    # get SCI value of final good group for semiramp
    data_final = data_masked[(ind_lastnz), range(data_masked.shape[1])]
    data_diff = data_final - data_zero # correctly does *NOT* have nans

    ind_lastnz = 0
    data_zero = 0

    # Use the readnoise and gain for good pixels only
    rn_sect_rav = rn_sect.flatten()[ good_pix ]
    rn_2_r = rn_sect_rav * rn_sect_rav

    gain_sect_r = gain_sect.flatten()[ good_pix ]

    # Calculate the sigma for nonzero gain values
    sigma_ir = data_final.copy() * 0.0
    numer_ir = data_final.copy() * 0.0
    data_final = 0

    # Calculate the SNR for pixels from the readnoise, the gain, and the
    # difference between the last and first reads for pixels where this results
    # in a positive SNR. Otherwise set the SNR to 0.
    sqrt_arg = rn_2_r + data_diff * gain_sect_r
    wh_pos = np.where((sqrt_arg >= 0.) & (gain_sect_r != 0.))
    numer_ir[wh_pos] = np.sqrt(rn_2_r[wh_pos] + \
                                data_diff[wh_pos] * gain_sect_r[wh_pos])
    sigma_ir[wh_pos] = numer_ir[wh_pos] / gain_sect_r[wh_pos]
    snr = data_diff * 0.
    snr[wh_pos] = data_diff[wh_pos] / sigma_ir[wh_pos]
    snr[snr < 0.] = 0.0

    del wh_pos

    gain_sect_r = 0
    numer_ir = 0
    data_diff = 0
    sigma_ir = 0

    power_wt_r = calc_power(snr)  # get the weighting exponent for this SNR

    # Make array of number of good groups, and exponents for each pixel
    num_nz = (data_masked != 0.).sum(0) # number of nonzero groups per pixel
    nrd_data_a = num_nz.copy()
    num_nz = 0

    nrd_prime = (nrd_data_a - 1) / 2.
    nrd_data_a = 0

    # Calculate inverse read noise^2 for use in weights
    invrdns2_r = 1./rn_2_r
    rn_sect = 0
    fnz = 0

    # Set optimal weights for each group of each pixel;
    #    for all pixels at once, loop over the groups
    wt_h = np.zeros(data_masked.shape, dtype=np.float32)

    for jj_rd in range(data_masked.shape[0]):
        wt_h[jj_rd, :] = \
            abs((abs(jj_rd - nrd_prime) / nrd_prime) ** power_wt_r) * invrdns2_r

    wt_h[np.isnan(wt_h)] = 0.
    wt_h[np.isinf(wt_h)] = 0.

    nrd_prime = 0
    power_wt_r = 0

    # For all pixels, 'roll' up the leading zeros such that the 0th group of
    #  each pixel is the lowest nonzero group for that pixel
    wh_m2d_f = np.logical_not(c_mask_2d[0, :]) # ramps with initial group False
    while (wh_m2d_f.sum() > 0):
        data_masked[:, wh_m2d_f] = np.roll(data_masked[:, wh_m2d_f], -1, axis=0)
        c_mask_2d[:, wh_m2d_f] = np.roll(c_mask_2d[:, wh_m2d_f], -1, axis=0)
        xvalues[:, wh_m2d_f] = np.roll(xvalues[:, wh_m2d_f], -1, axis=0)
        wh_m2d_f = np.logical_not(c_mask_2d[0, :])

    # Create weighted sums for Poisson noise and read noise
    nreads_wtd = (wt_h * c_mask_2d).sum(axis=0)  # using optimal weights

    sumx = (xvalues * wt_h).sum(axis=0)
    sumxx = (xvalues**2 * wt_h).sum(axis=0)

    c_data_masked = data_masked.copy()
    c_data_masked[np.isnan(c_data_masked)] = 0.
    sumy = (np.reshape((c_data_masked * wt_h).sum(axis=0), sumx.shape))
    sumxy = (xvalues * wt_h * np.reshape(c_data_masked, xvalues.shape)).sum(axis=0)

    return sumx, sumxx, sumxy, sumy, nreads_wtd, xvalues


def log_stats(c_rates):
    """
    Short Summary
    -------------
    Optionally log statistics of detected cosmic rays

    Parameters
    ----------
    c_rates: float, 2D array
       weighted count rate

    Returns
    -------
    None
    """
    wh_c_0 = np.where(c_rates == 0.) # insuff data or no signal

    log.debug('The number of pixels having insufficient data')
    log.debug('due to excessive CRs or saturation %d:', len(wh_c_0[0]))
    log.debug('Count rates - min, mean, max, std: %f, %f, %f, %f'
             % (c_rates.min(), c_rates.mean(), c_rates.max(), c_rates.std()))

