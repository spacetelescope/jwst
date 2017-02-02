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
#  In this module, comments on the 'first read','second read', etc are 1-based.

from __future__ import division
import time
import numpy as np
import logging

from .. import datamodels
from ..datamodels import dqflags

from . import gls_fit           # used only if algorithm is "GLS"
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

BUFSIZE = 1024 * 30000  # 30Mb cache size for data section
MIN_ERR = 1.E-10 # minimum allowed value for error cube
MIN_LEN = 1 # MIN_LEN+1 is min length of intervals to be fit by least squares

# Replace zero or negative variances with this.
LARGE_VARIANCE = 1.e8


def ramp_fit(model, buffsize, save_opt, readnoise_model, gain_model,
              algorithm, weighting):
    """
    Extended Summary
    ----------------
    Calculate the count rate for each pixel in all data cube sections and all
    integrations, equal to the slope for all sections (intervals between
    cosmic rays) of the pixel's ramp divided by the effective integration time.
    If the weighting parameter is set to 'optim', the optimal weighting from a
    pere by Fixsen (ref. TBA) will be used in the fitting; otherwise the
    fitting will be unweighted.

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
        'unweighted' specifies that no weighting should be used (default)
        'optimal' specifies that optimal weighting should be used

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
        new_model, int_model, opt_model = ols_ramp_fit(model,
                                buffsize, save_opt,
                                readnoise_model, gain_model, weighting)
        gls_opt_model = None


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
        'unweighted' specifies that no weighting should be used (default)
        'optimal' specifies that optimal weighting should be used

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

    # get needed sizes and shapes
    nreads, npix, imshape, cubeshape, n_int, instrume, frame_time, ngroups = \
        utils.get_dataset_info(model)

    # Save original shapes for writing to log file, as these may change for MIRI
    orig_nreads = nreads
    orig_cubeshape = cubeshape

    # For MIRI datasets having >1 reads, if all final reads are flagged
    # as DO_NOT_USE, resize the input model arrays to exclude the final group.
    if (instrume == 'MIRI' and nreads > 1):
        last_gdq = model.groupdq[:,-1,:,:]
        gdq_shape = last_gdq.shape

        if np.array_equal(np.full(gdq_shape, dqflags.group['DO_NOT_USE']), \
                      last_gdq):
            model.data = model.data[:,:-1,:,:]
            model.err = model.err[:,:-1,:,:]
            model.groupdq = model.groupdq[:,:-1,:,:]
            nreads -= 1
            ngroups -= 1
            cubeshape = (nreads,)+imshape
            log.info('MIRI dataset has all final reads flagged as DO_NOT_USE.')

    if (ngroups == 1):
        log.warn('Dataset has NGROUPS=1, so count rates for each integration')
        log.warn('will be calculated as the value of that 1 group divided by')
        log.warn('the group exposure time.')

    slopes = np.zeros(imshape, dtype=np.float32)

    slope_wtd = np.zeros(imshape, dtype=np.float64)
    m_sum_2d = np.zeros(imshape, dtype=np.float64)
    var_sum_2d = np.zeros(imshape, dtype=np.float64)

    # For multiple-integration datasets, will output integration-specific
    #    results to separate file named <basename> + '_integ.fits'
    slope_int, err_int, dq_int, m_by_var_int, inv_var_int = \
             utils.alloc_int(n_int, imshape)

    # Get GROUP DQ and ERR arrays from input file
    gdq_cube = model.groupdq
    err_cube = model.err

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
    readnoise_2d, gain_2d = utils.get_ref_subs(model, readnoise_model, gain_model)

    if save_opt:
            # get max number of segments fit in all integrations
        max_seg = calc_num_seg(gdq_cube, imshape, n_int, nreads)

        # create object to hold optional, segment-specific results
        opt_res = utils.OptRes(n_int, imshape, max_seg, nreads)
    else:
        max_seg = 1 # needed for calc_slope()
        opt_res = None

    # loop over data integrations
    for num_int in range(0, n_int):

        # loop over data sections
        for rlo in range(0, cubeshape[1], nrows):
            rhi = rlo + nrows

            if rhi > cubeshape[1]:
                rhi = cubeshape[1]

            data_sect = model.get_section('data')[num_int, :, rlo:rhi, :]

            # first frame section for 1st read of current integration
            ff_sect = model.get_section('data')[num_int,
                                                0, rlo:rhi, :].astype(np.float32)
            # get appropriate sections
            gdq_sect = gdq_cube[num_int, :, rlo:rhi, :]
            rn_sect = readnoise_2d[rlo:rhi, :]
            gain_sect = gain_2d[rlo:rhi, :]

            t_err_cube, t_dq_cube, m_by_var, inv_var, opt_res = \
                 calc_slope(data_sect, gdq_sect, frame_time, opt_res, \
                             rn_sect, gain_sect, max_seg, ngroups, weighting)

            err_cube[num_int, :, rlo:rhi, :] += t_err_cube
            gdq_cube[num_int, :, rlo:rhi, :] = t_dq_cube

            # Compress 4D->2D dq arrays for saturated and jump-detected pixels
            pixeldq_sect = pixeldq[rlo:rhi, :].copy()
            dq_int[num_int, rlo:rhi, :] = \
                  dq_compress_sect(t_dq_cube, pixeldq_sect).copy()

            sect_shape = data_sect.shape[-2:]
            m_sum_2d[rlo:rhi, :] += m_by_var.reshape(sect_shape)
            var_sum_2d[rlo:rhi, :] += inv_var.reshape(sect_shape)

            if save_opt: # collect optional results for output
                opt_res.reshape_res(num_int, rlo, rhi, sect_shape, ff_sect)

            # Calculate difference between each slice and the previous slice
            #   as approximation to cosmic ray amplitude for those pixels
            #   having their DQ set for cosmic rays
                data_diff = data_sect - utils.shift_z(data_sect, -1)
                dq_cr = np.bitwise_and(dqflags.group['JUMP_DET'], gdq_sect)

                opt_res.cr_mag_seg[num_int, :, rlo:rhi, :] = data_diff * (dq_cr != 0)

            m_by_var_int[num_int, rlo:rhi, :] = m_by_var.reshape(sect_shape)
            inv_var_int[num_int, rlo:rhi, :] = inv_var.reshape(sect_shape)

        slope_int[num_int, :, :] = \
                 utils.calc_slope_int(slope_int, m_by_var_int, inv_var_int,
                                       num_int)

        err_int[num_int, :, :] = err_cube[0, 0, :, :] # change eventually

        if save_opt: # collect optional pedestal results for output
            opt_res.ped_int[num_int, :, :] = \
                   utils.calc_pedestal(num_int, slope_int, opt_res.firstf_int,
                                       gdq_cube)

    wh_non_zero = (var_sum_2d != 0.0)

    slope_wtd[wh_non_zero] = (m_sum_2d[wh_non_zero] /
                                var_sum_2d[wh_non_zero])

    slopes = slope_wtd.reshape(imshape)


    # Calculate effective integration time (once EFFINTIM has been populated
    #   and accessible, will use that instead)
    effintim = utils.get_effintim(model)

    if save_opt: # collect optional results for output
        ### opt_res.print_full() # diagnostic; uncomment for small datasets
        # Shrink cosmic ray magnitude array; attach to optional results object
        opt_res.shrink_crmag(n_int, gdq_cube, imshape, nreads)
        opt_model = opt_res.output_optional(model, effintim)
    else:
        opt_model = None

    # Divide slopes by total (summed over all integrations) effective
    #   integration time to give count rates.
    c_rates = slopes / effintim

    # Values in err_cube have been 1./variance, so take reciprocal
    #    of non-zero values to write
    err_cube[err_cube <= 0.] = MIN_ERR  # for pixels having no signal
    err_cube = 1. / err_cube # has no zero values

    # Compress all integration's dq arrays to create 2D PIXELDDQ array for
    #   primary output
    final_pixeldq = dq_compress_final(dq_int, n_int)

    if n_int > 1:
        int_model = utils.output_integ(model, slope_int, err_int, dq_int,
                                        effintim)
    else:
        int_model = None

    tstop = time.time()

    log_stats(c_rates)

    log.debug('Instrument: %s' % (instrume))
    log.debug('Number of pixels in 2D array: %d' % (npix))
    log.debug('Shape of 2D image: (%d, %d)' % (imshape))
    log.debug('Shape of data cube: (%d, %d, %d)' % (orig_cubeshape))
    log.debug('Buffer size (bytes): %d' % (buffsize))
    log.debug('Number of rows per buffer: %d' % (nrows))
    log.info('Number of groups per integration: %d' % (orig_nreads))
    log.info('Number of integrations: %d' % (n_int))
    log.debug('The execution time in seconds: %f' % (tstop - tstart))

    # Create new model...
    new_model = datamodels.ImageModel(data=c_rates.astype(np.float32),
                                   dq=final_pixeldq.astype(np.int32),
                                   err=err_cube[0, 0].copy())

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
        slopes = np.zeros(imshape, dtype=np.float64)
        sum_weight = np.zeros(imshape, dtype=np.float64)

    # For multiple-integration datasets, will output integration-specific
    # results to separate file named <basename> + '_integ.fits'.
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
        i_max_num_cr = utils.get_max_num_cr(gdq_cube[num_int, :, :, :], jump_flag)
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

    # These parameters will be compared with both the readnoise and
    # gain models.
    xstart = model.meta.subarray.xstart
    xsize = model.meta.subarray.xsize
    ystart = model.meta.subarray.ystart
    ysize = model.meta.subarray.ysize

    # Get readnoise array.
    if (readnoise_model.meta.subarray.xstart == xstart and
        readnoise_model.meta.subarray.xsize == xsize and
        readnoise_model.meta.subarray.ystart == ystart and
        readnoise_model.meta.subarray.ysize == ysize):

        log.debug('Readnoise subarray matches science data')
        readnoise_2d = readnoise_model.data
    else:
        log.debug('Extracting readnoise subarray to match science data')
        xstop = xstart + xsize - 1
        ystop = ystart + ysize - 1
        readnoise_2d = readnoise_model.data[ystart - 1:ystop, xstart - 1:xstop]

    # Get gain array.
    if (gain_model.meta.subarray.xstart == xstart and
        gain_model.meta.subarray.xsize == xsize and
        gain_model.meta.subarray.ystart == ystart and
        gain_model.meta.subarray.ysize == ysize):

        log.debug('Gain subarray matches science data')
        gain_2d = gain_model.data
    else:
        log.debug('Extracting gain subarray to match science data')
        xstop = xstart + xsize - 1
        ystop = ystart + ysize - 1
        gain_2d = gain_model.data[ystart - 1:ystop, xstart - 1:xstop]

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
    log.info('Number of reads per integration: %d' % nreads)
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
                                      dq=final_pixeldq,
                                      err=gls_err.astype(np.float32))
    else:
        new_model = datamodels.ImageModel(data=slope_int[0],
                                      dq=final_pixeldq,
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
    Get ramp locations where the data has been flagged as saturated in the
    4D GROUPDQ array for the current data section, find the corresponding image
    locations, and set the SATURATED flag in those locations in the
    PIXELDQ array. Similarly, get the ramp locations where the data has been
    flagged as a jump detection in the 4D GROUPDQ array, find the corresponding
    image locations, and set the COSMIC_BEFORE flag in those locations in the
    PIXELDQ array

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


def calc_slope(data_sect, gdq_sect, frame_time, opt_res, rn_sect, gain_sect,
                max_seg, ngroups, weighting):
    """
    Short Summary
    -------------
    Calculate the count rate for each pixel in the data cube section
    for the current integration, equal to the weighted slope for all
    segments (intervals between cosmic rays) of the pixel's ramp.

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

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    max_seg: int
        maximum number of segments that will be fit within an
        integration, calculated over all pixels and all integrations

    ngroups: int
        number of groups per integration

    weighting: string
        'unweighted' specifies that no weighting should be used (default)
        'optimal' specifies that optimal weighting should be used

    Returns
    -------
    err_sect: float, 3D array
        fitting error estimate for pixels in section

    gdq_sect: int, 3D array
        data quality flags for pixels in section

    m_by_var: float, 1D array
        values of slope/variance for good pixels

    inv_var: float, 1D array
        values of 1/variance for good pixels

    opt_res: OptRes object
        contains quantities related to fitting for optional output

    """
    nreads, asize2, asize1 = data_sect.shape
    npix = asize2 * asize1  # number of pixels in section of 2D array
    imshape = data_sect.shape[-2:]
    cubeshape = (nreads,) + imshape  # cube section shape

    all_pix = np.arange(npix)

    arange_nreads_col = np.arange(nreads)[:, np.newaxis]

    start = np.zeros(npix, dtype=np.int32) # lowest channel in fit

    # Highest channel in fit initialized to last read
    end = np.zeros(npix, dtype=np.int32) + (nreads - 1)

    pixel_done = (end < 0) # False until processing is done

    # Weighted average for slopes, equal to sum(m/v)/sum(1/v)
    #    for the pixels in which the variance v is nonzero
    inv_var = np.zeros(npix, dtype=np.float64)
    m_by_var = np.zeros(npix, dtype=np.float64)

    num_seg = np.zeros(npix, dtype=np.int32)

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

    wh_f = np.where(mask_2d == False)
    these_p = wh_f[1] # coordinates of pixels flagged as False
    these_r = wh_f[0] # reads of pixels flagged as False

    total_mask_sum = mask_2d.sum(axis=0) # initial number of good reads along ramp

    # Populate end_st to contain the set of end points for each pixel.
    # Populate end_heads to initially include every pixel that is either
    # saturated or contains a cosmic ray. Skips the duplicated final read
    # for saturated pixels. Saturated pixels resulting in a contiguous set
    # of intervals of length 1 will later be flagged as too short
    # to fit well.
    for ii in range(len(these_p)):
        if (these_r[ii] != (nreads - 1)):
            end_st[end_heads[these_p[ii]], these_p[ii]] = these_r[ii]
            end_heads[these_p[ii]] += 1

    # Sort and reverse array to handle the order that saturated pixels
    # were added
    end_st.sort(axis=0)
    end_st = end_st[::-1]

    # Reformat to designate the number of endpoints per pixel; compress
    # to specify number of reads per pixel
    end_heads = (end_st > 0).sum(axis=0)

    # Create object to hold optional results
    if (opt_res is not None):
        opt_res.init_2d(npix, max_seg)

    # LS fit until 'nreads' iterations or all pixels in
    #    section have been processed
    for iter_num in range(nreads):
        if pixel_done.all():
            break

        # frames >= start and <= end_st will be included in fit
        mask_2d = ((arange_nreads_col >= start) &
                   (arange_nreads_col <
                    (end_st[end_heads[all_pix] - 1, all_pix] + 1)))

        mask_2d[gdq_sect_r != 0] = False # exclude bad group dq values

        # for all pixels, update arrays, summing slope and variance
        fit_next_segment(start, end_st, end_heads, pixel_done, data_sect,
                          mask_2d, inv_var, m_by_var, num_seg, opt_res,
                          rn_sect, gain_sect, ngroups, weighting, total_mask_sum)

    arange_nreads_col = 0
    all_pix = 0

    err_sect = np.zeros(cubeshape, dtype=np.float32)

    # For now, making all error array slices within an integration and
    #  section identical. Update later when use of error array has been decided
    for ii in range(cubeshape[0]):
        err_sect[ii, :, :] = err_2d_array

    return err_sect, gdq_sect, m_by_var, inv_var, opt_res


def fit_next_segment(start, end_st, end_heads, pixel_done, data_sect, mask_2d,
                      inv_var, m_by_var, num_seg, opt_res, rn_sect, gain_sect,
                      ngroups, weighting, total_mask_sum):
    """
    Extended Summary
    ----------------
    Call routine to LS fit masked data for a single segment for all
    pixels in data section. Then categorize each pixel's fitting
    interval based on interval length, and whether the interval is at
    the end of the array.  Update the start array, the end stack array,
    the end_heads array which contains the number of endpoints. For
    pixels in which the fitting intervals are long enough, the
    resulting slope and variance are added to the appropriate stack
    arrays.  The first channel to fit in a segment is either the first
    read in the ramp, or a read in which a cosmic ray has been flagged.

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

    mask_2d: int, 2D array
        delineates which channels to fit for each pixel

    inv_var: float, 1D array
        values of 1/variance for good pixels

    m_by_var: float, 1D array
        values of slope/variance for good pixels

    num_seg: int, 1D array
        numbers of segments for good pixels

    opt_res: OptRes object
        optional fitting results to output

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    ngroups: int
        number of groups per integration

    weighting: string
        'unweighted' specifies that no weighting should be used (default)
        'optimal' specifies that optimal weighting should be used

    total_mask_sum: int, 1D array
        initial number of good reads along ramp, used to distinguish between
        cosmic rays and saturated pixels

    Returns
    -------
    None

    """
    nreads, asize2, asize1 = data_sect.shape # Note: nreads is a scalar here
    all_pix = np.arange(asize2 * asize1)

    slope, intercept, variance, sig_intercept, sig_slope = \
        fit_lines(data_sect, mask_2d, end_heads, end_st, start, rn_sect,
                   gain_sect, ngroups, weighting)

    end_locs = end_st[end_heads[all_pix] - 1, all_pix]
    l_interval = end_locs - start # fitting interval length

    wh_done = (start == -1) # done pixels
    l_interval[wh_done] = 0  # set interval lengths for done pixels to 0

    # CASE 0 - full-length ramp has 1 good read, and saturation on 2nd read
    #    or all reads flagged on 2nd to end, so will use single good read
    #    data as the slope
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of end to 0
    #    - set pixel_done to True to designate all fitting done
    wh_sat2 = np.where((total_mask_sum == 1) & (l_interval == 1) &
              (end_locs != nreads - 1) & (mask_2d[:, :].sum(axis=0) == 1))

    if(len(wh_sat2[0]) > 0):
        these_pix = wh_sat2[0]

        start[these_pix] = -1
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] = 0
        pixel_done[these_pix] = True # all processing for pixel is completed
        inv_var[these_pix] += 1.0 / variance[these_pix]
        m_by_var[these_pix] += slope[these_pix] / variance[these_pix]

        if (opt_res is not None):
            # Append results to arrays
            opt_res.append_arr(num_seg, these_pix, intercept, slope,\
                                    sig_intercept, sig_slope, inv_var)

        num_seg[these_pix] += 1


    # CASE 0.5 interval is the first segment in a ramp having a good 1st read,
    #     a cosmic ray in the 2nd read, and at least one good read later.
    #    - set start to 1 beyond end of current interval
    #    - remove current end from end stack
    #    - decrement number of ends

    wh_sat2 = np.where((total_mask_sum > 1) & (l_interval == 1) &
               (end_locs != nreads - 1) & (mask_2d[:, :].sum(axis=0) == 1))

    if(len(wh_sat2[0]) > 0):
        these_pix = wh_sat2[0]

        start[these_pix] = end_locs[these_pix]
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] -= 1
        wh_neg = (end_heads < 0.)
        end_heads[wh_neg] = 0.


    # CASE 1 - interval too short to fit well, not at array end,
    #    (but not saturated on either 1st or 2nd read (this latter is case 0)
    #    - set start to 1 beyond end of current interval
    #    - remove current end from end stack
    #    - decrement number of ends
    wh_check = np.where(((l_interval <= MIN_LEN) & (end_locs != nreads - 1)) &
                       ~((l_interval == 1) & (mask_2d[:, :].sum(axis=0) == 1)))
    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        start[these_pix] = end_locs[these_pix]
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] -= 1
        wh_neg = (end_heads < 0.)
        end_heads[wh_neg] = 0.

    # CASE 2 - interval too short to fit well, at end of array, NRGOUPS>1,
    #    but exclude NGROUPS==2 datasets as they are covered in case 6.
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - set pixel_done to True to designate all fitting done
    wh_check = np.where((l_interval <= MIN_LEN) & (end_locs == nreads - 1) &\
                        (nreads > 1) & (ngroups != 2))
    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        start[these_pix] = -1
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] = 0
        pixel_done[these_pix] = True

    # CASE 3 - interval long enough, not at end of array
    #    - remove current end from end stack
    #    - decrement number of ends
    #    - add slopes and variances to running sums
    wh_check = np.where((l_interval > MIN_LEN) & (end_locs != nreads - 1))
    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        start[these_pix] = end_locs[these_pix]
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] -= 1
        wh_neg = (end_heads < 0.)
        end_heads[wh_neg] = 0.

        g_pix = these_pix[variance[these_pix] > 0.] # good pixels

        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]
            m_by_var[g_pix] += slope[g_pix] / variance[g_pix]

            if (opt_res is not None):
                # Append results to arrays
                opt_res.append_arr(num_seg, g_pix, intercept, slope, \
                                    sig_intercept, sig_slope, inv_var)

            num_seg[g_pix] += 1

    # CASE 4 - interval long enough, at end of array
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - add slopes and variances to running sums
    wh_check = np.where((l_interval > MIN_LEN) & (end_locs == nreads - 1))
    if(len(wh_check[0]) > 0):
        these_pix = wh_check[0]
        start[these_pix] = -1   # all processing for this pixel is completed
        end_st[end_heads[these_pix] - 1, these_pix] = 0
        end_heads[these_pix] = 0
        pixel_done[these_pix] = True # all processing for pixel is completed

        g_pix = these_pix[variance[these_pix] > 0.] # good pixels
        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]
            m_by_var[g_pix] += slope[g_pix] / variance[g_pix]

            if (opt_res is not None):
                # Append results to arrays
                opt_res.append_arr(num_seg, g_pix, intercept, slope,\
                                    sig_intercept, sig_slope, inv_var)

            num_seg[g_pix] += 1

    # CASE 5 - NGROUPS=1 or NGROUPS=2; so special fitting is done for all pixels,
    #    and all intervals are at the end of the array. Otherwise, this Case is
    #    identical to Case 2.
    #    - set start to -1 to designate all fitting done
    #    - remove current end from end stack
    #    - set number of ends to 0
    #    - set pixel_done to True to designate all fitting done
    if (ngroups == 1 or ngroups == 2):
        start[all_pix] = -1
        end_st[end_heads[all_pix] - 1, all_pix] = 0
        end_heads[all_pix] = 0
        pixel_done[all_pix] = True

        g_pix = all_pix[variance[all_pix] > 0.]

        if (len(g_pix) > 0):
            inv_var[g_pix] += 1.0 / variance[g_pix]
            m_by_var[g_pix] += slope[g_pix] / variance[g_pix]

            if (opt_res is not None): # multi-integration, NGROUPS=1 datasets
                # Append results to arrays
                opt_res.append_arr(num_seg, g_pix, intercept, slope,\
                                    sig_intercept, sig_slope, inv_var)

            num_seg[g_pix] = 1

    return None


def fit_lines(data, mask_2d, end_heads, end_st, start, rn_sect, gain_sect,
                ngroups, weighting):
    """
    Extended Summary
    ----------------
    Do linear least squares fit to data cube in this integration.  This
    function will later be generalized to accept input data of arbitrary
    dimensionality. In addition to applying the mask due to identified
    cosmic rays, the data is also masked to exclude intervals that are
    too short to fit well. The first channel to fit in a segment is either
    the first read in the ramp, or a read in which a cosmic ray has been
    flagged.

    Parameters
    ----------
    data: float
       array of values for current data section

    mask_2d: boolean, 2D array
       delineates which channels to fit for each pixel

    end_heads: int, 1D array
       number of endpoints for each pixel

    end_st: int, 2D array
       stack array of endpoints

    start: int, 1D array
        lowest channel in fit

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    ngroups: int
        number of groups per integration

    weighting: string
        'unweighted' specifies that no weighting should be used (default)
        'optimal' specifies that optimal weighting should be used

    Returns
    -------
    full_slope: float, 1D array
       weighted slope for current iteration's pixels for data section

    full_intercept: float, 1D array
       y-intercepts from fit for data section

    full_variance: float, 1D array
       variance of residuals for fit for data section

    full_sig_intercept: float, 1D array
       sigma of y-intercepts from fit for data section

    full_sig_slope: float, 1D array
       sigma of slopes from fit for data section

    """
    # verify that incoming data is either 2 or 3-dimensional
    try:
        assert (data.ndim == 2 or data.ndim == 3)
    except AssertionError:
        log.error('FATAL ERROR: Data input to fit_lines must be 2 or 3 dimensions')

    # To ensure that the first channel to be fit is the cosmic-ray-affected
    #   read, the channel previous to each channel masked as good is
    #   also masked as good. This is only for the local purpose of setting
    #   the first channel, and will not propagate beyond this current function
    #   call.
    wh_mask_2d = np.where(mask_2d)
    mask_2d[np.maximum(wh_mask_2d[0] - 1, 0), wh_mask_2d[1]] = True

    nreads_1d = mask_2d.astype(np.int16).sum(axis=0) # num of reads/pixel unmasked

    npix = mask_2d.shape[1]

    full_slope = np.zeros(npix, dtype=np.float64)
    full_variance = np.zeros(npix, dtype=np.float64) + MIN_ERR
    full_intercept = np.zeros(npix, dtype=np.float64)
    full_sig_intercept = np.zeros(npix, dtype=np.float64) + MIN_ERR
    full_sig_slope = np.zeros(npix, dtype=np.float64) + MIN_ERR

    # For full-length ramps in the current data section in which only the
    # 1st read has good data (these are ramps in which the 2nd read is
    # saturated) set the slope equal to that good data value, and set the
    # intercept to 0.  This is done before the arrays are compressed.
    wh_sat1 = np.where((mask_2d[:, :].sum(axis=0) == 1) & (mask_2d[0, :] == True)
                       & (end_heads == data.shape[0] - 1))
    if (len(wh_sat1[0]) > 0):
        data0_slice = data[0, :, :].reshape(npix)
        full_slope[wh_sat1] = data0_slice[wh_sat1]
        full_variance[wh_sat1] = MIN_ERR
        full_sig_slope[wh_sat1] = MIN_ERR
        full_intercept[wh_sat1] = 0.
        full_sig_intercept[wh_sat1] = 0.

    all_pix = np.arange(npix)
    end_locs = end_st[end_heads[all_pix] - 1, all_pix]
    l_interval = end_locs - start # fitting interval length

    mask_2d[:, l_interval < MIN_LEN] = False

    wh_pix_to_use = np.where(mask_2d.sum(axis=0) > MIN_LEN)
    good_pix = wh_pix_to_use[0]

    # reshape data_masked and eliminate pixels with too-short intervals
    data_masked = data * np.reshape(mask_2d, data.shape)
    data_masked = np.reshape(data_masked, (data_masked.shape[0], npix))
    data_masked = data_masked[:, good_pix]

    xvalues = np.arange(data_masked.shape[0])[:, np.newaxis] * mask_2d  # all
    xvalues = xvalues[:, good_pix]  # set to those pixels to be used

    mask_2d = mask_2d[:, good_pix]
    nreads_1d = nreads_1d[good_pix]

    if weighting.lower() == 'optimal': # do the fits using optimal weighting
        # get sums from optimal weighting
        sumx, sumxx, sumxy, sumy, nreads_wtd, xvalues =\
           calc_opt_sums(rn_sect, gain_sect, data_masked, mask_2d, \
                          xvalues)

        slope, intercept, sig_slope, sig_intercept =\
               calc_opt_fit(nreads_wtd, sumxx, sumx, sumxy, sumy)

        denominator = nreads_wtd * sumxx - sumx**2

    elif weighting.lower() == 'unweighted': # do the fits using unweighted weighting
        # get sums from unweighted weighting
        sumx, sumxx, sumxy, sumy =\
              calc_unwtd_sums(data_masked, xvalues)

        denominator = nreads_1d * sumxx - sumx**2

        slope, intercept, sig_slope, sig_intercept, line_fit =\
               calc_unwtd_fit(xvalues, nreads_1d, sumxx, sumx, sumxy, sumy)

    else: # unsupported weighting type specified
        log.error('FATAL ERROR: unsupported weighting type specified.')

    line_fit = 0
    variance = nreads_1d / denominator
    denominator = 0

    # check to prevent NaN propagation
    variance = correct_noiseless(variance, nreads_1d, data, good_pix, rn_sect)

    full_slope[good_pix] = slope
    full_variance[good_pix] = variance
    full_intercept[good_pix] = intercept
    full_sig_intercept[good_pix] = sig_intercept
    full_sig_slope[good_pix] = sig_slope

    # include pixels with >1 good reads
    mask_2d = mask_2d.compress(mask_2d.sum(axis=0) > 1., axis=1)

    if (ngroups == 1): # process 1 group/integration dataset
        full_slope, full_intercept, full_variance, full_sig_intercept, \
        full_sig_slope = fit_1_group(full_slope, full_intercept, \
        full_variance, full_sig_intercept, full_sig_slope, npix, data, mask_2d)

    if (ngroups == 2): # process 2 group/integration dataset
        full_slope, full_intercept, full_variance, full_sig_intercept, \
        full_sig_slope = fit_2_group(full_slope, full_intercept, \
        full_variance, full_sig_intercept, full_sig_slope, npix, data, mask_2d)

    return full_slope, full_intercept, full_variance,  \
           full_sig_intercept, full_sig_slope


def calc_unwtd_fit(xvalues, nreads_1d, sumxx, sumx, sumxy, sumy):
    """
    Extended Summary
    ----------------
    Do linear least squares fit to data cube in this integration, using
    unweighted fits to the segments.

    Parameters
    ----------
    xvalues: int, 1D array
        indices of valid pixel values for all reads

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
    Do linear least squares fit to data cube in this integration, using
    optimally weighted fits to the segments.  The weighting uses the formulation
    by Fixsen (Fixsen et al, PASP, 112, 1350).

    Parameters
    ----------
    nreads_wtd: float, 1D array
        sum of product of data and weight

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
    sig_slope = (nreads_wtd / denominator)**0.5

    # Set to 0 those values that are NaN or Inf just in case these were not
    # properly handled by dq
    slope[np.isnan(slope)] = 0.
    intercept[np.isnan(intercept)] = 0.
    sig_slope[np.isnan(sig_slope)] = 0.
    sig_intercept[np.isnan(sig_intercept)] = 0.

    return slope, intercept, sig_slope, sig_intercept



def fit_1_group(full_slope, full_intercept, full_variance, full_sig_intercept,
                 full_sig_slope, npix, data, mask_2d):
    """
    Extended Summary
    ----------------
    This function sets the fitting arrays for datasets having only 1 group
    per integration.

    Parameters
    ----------
    full_slope: float, 1D array
        weighted slope for current iteration's pixels for data section

    full_intercept: float, 1D array
        y-intercepts from fit for data section

    full_variance: float, 1D array
        variance of residuals for fit for data section

    full_sig_intercept: float, 1D array
        sigma of y-intercepts from fit for data section

    full_sig_slope: float, 1D array
        sigma of slopes from fit for data section

    npix: int
        number of pixels in 2d array

    data: float
        array of values for current data section

    mask_2d: int, 2D array
        delineates which channels to fit for each pixel

    Returns
    -------
    full_slope: float, 1D array
        weighted slope for current iteration's pixels for data section

    full_intercept: float, 1D array
        y-intercepts from fit for data section

    full_variance: float, 1D array
        variance of residuals for fit for data section

    full_sig_intercept: float, 1D array
        sigma of y-intercepts from fit for data section

    full_sig_slope: float, 1D array
        sigma of slopes from fit for data section
    """
    # For pixels not saturated, recalculate the slope as the value of the SCI
    # data in that group, which will later be divided by the group exposure
    # time to give the count rate. Recalculate other fit quantities to be
    # benign.
    full_slope = data[0, :, :].reshape(npix)
    full_variance = np.zeros(npix, dtype=np.float64) + MIN_ERR
    full_sig_slope = full_slope * 0. + MIN_ERR
    full_intercept = full_slope * 0.
    full_sig_intercept = full_slope * 0. + MIN_ERR

    # For saturated pixels, overwrite fit values with benign values.
    wh_sat0 = np.where(mask_2d[0, :] == False)
    if (len(wh_sat0[0]) > 0):
        sat_pix = wh_sat0[0]
        full_slope[sat_pix] = 0.
        full_variance[sat_pix] = MIN_ERR
        full_sig_slope[sat_pix] = MIN_ERR
        full_intercept[sat_pix] = 0.
        full_sig_intercept[sat_pix] = MIN_ERR

    return full_slope, full_intercept, full_variance, full_sig_intercept, full_sig_slope


def fit_2_group(full_slope, full_intercept, full_variance, full_sig_intercept,
                 full_sig_slope, npix, data, mask_2d):
    """
    Extended Summary
    ----------------
    This function sets the fitting arrays for datasets having only 2 groups
    per integration.

    Parameters
    ----------
    full_slope: float, 1D array
        weighted slope for current iteration's pixels for data section

    full_intercept: float, 1D array
        y-intercepts from fit for data section

    full_variance: float, 1D array
        variance of residuals for fit for data section

    full_sig_intercept: float, 1D array
        sigma of y-intercepts from fit for data section

    full_sig_slope: float, 1D array
        sigma of slopes from fit for data section

    npix: int
        number of pixels in 2d array

    data: float
        array of values for current data section

    mask_2d: int, 2D array
        delineates which channels to fit for each pixel

    Returns
    -------
    full_slope: float, 1D array
        weighted slope for current iteration's pixels for data section

    full_intercept: float, 1D array
        y-intercepts from fit for data section

    full_variance: float, 1D array
        variance of residuals for fit for data section

    full_sig_intercept: float, 1D array
        sigma of y-intercepts from fit for data section

    full_sig_slope: float, 1D array
        sigma of slopes from fit for data section
    """
    # For pixels saturated on the first read, overwrite fit values with
    # benign values.
    wh_sat0 = np.where(mask_2d[0, :] == False)
    if (len(wh_sat0[0]) > 0):
        sat_pix = wh_sat0[0]
        full_slope[sat_pix] = 0.
        full_variance[sat_pix] = MIN_ERR
        full_sig_slope[sat_pix] = MIN_ERR
        full_intercept[sat_pix] = 0.
        full_sig_intercept[sat_pix] = 0.

    # For pixels saturated on the second read, recalculate the slope as
    # the value of the SCI data in the first read, which will later be
    # divided by the group exposure time to give the count rate, and
    # recalculate the other fit quantities to be benign. Note: these pixels
    # will already have been handled earlier (for intervals of arbitrary
    # length) in this function, but are being included here to explicitly
    # cover all possibilities for pixels in datasets with ngroups=2. Will
    # later consider refactoring.
    wh_sat1 = np.where((mask_2d[:, :].sum(axis=0) == 1) & (mask_2d[0, :] == True))
    if (len(wh_sat1[0]) > 0):
        data0_slice = data[0, :, :].reshape(npix)
        full_slope[wh_sat1] = data0_slice[wh_sat1]
        full_variance[wh_sat1] = MIN_ERR
        full_sig_slope[wh_sat1] = MIN_ERR
        full_intercept[wh_sat1] = 0.
        full_sig_intercept[wh_sat1] = 0.

    # For pixels with no saturated values, recalculate the slope as the
    # difference between the values of the second and first reads, which will
    # later be divided by the group exposure time to give the count rate, and
    # recalculate other fit quantities to be benign.
    wh_sat_no = np.where(mask_2d[:, :].sum(axis=0) == 2)
    if (len(wh_sat_no[0]) > 0):
        data0_slice = data[0, :, :].reshape(npix)
        data1_slice = data[1, :, :].reshape(npix)
        full_slope[wh_sat_no] = data1_slice[wh_sat_no] - data0_slice[wh_sat_no]
        ## full_variance - already been correctly calculated
        full_sig_slope[wh_sat_no] = MIN_ERR
        full_intercept[wh_sat_no] = data0_slice[wh_sat_no] -\
            data1_slice[wh_sat_no] # by geometry
        full_sig_intercept[wh_sat_no] = MIN_ERR

    return full_slope, full_intercept, full_variance, full_sig_intercept, \
           full_sig_slope


def correct_noiseless(variance, nreads, data, good_pix, rn_sect):
    """
    Short Summary
    -------------
    Semi-ramps that are either noiseless, or have no signal, or contain only
    2 reads will have variance = 0.  All such semi-ramps will have their
    variance recalculated here, equal to the poisson noise of the ramp added in
    quadrature to the read noise, ensuring that all variance values are
    positive.

    Parameters
    ----------
    variance: float, 1D array
        inverse weighting of semi-ramps

    nreads: float, 1D array
        number of reads in an integration

    data: float, 3D array
        SCI data for all pixels in semi-ramp in data section

    good_pix: int, 1D array
        indices of pixels having valid data for all reads

    rn_sect: float, 2D array
        read noise values for all pixels in data section

    Returns
    -------
    variance: float, 1D array
        variances of all good pixels, including the recalculated variances of
        linear fits of semi-ramps having incoming variances that are 0.

    """
    wh_varnan = np.isnan(variance)
    wh_numer0 = nreads == 0.0

    # Return immediately if there are no good pixels, or if there are no
    #   variance=0 pixels
    if (len(good_pix) == 0):
        return variance
    if ((wh_numer0.sum() + wh_varnan.sum()) == 0):
        return variance

    # There is at least 1 pixel having a noiseless semi-ramp, so calculate the
    #   poisson error of those semi-ramps, using the last read minus the first
    #   read
    data_diff = (data[-1, :, :] - data[0, :, :]).ravel()[good_pix]
    poiss_2d = np.reshape(data_diff, variance.shape)

    rn_2d = rn_sect.ravel()[good_pix]
    new_var = np.sqrt(poiss_2d + rn_2d**2.)

    # Overwrite pixels having noiseless semi-ramps with revised values
    variance[wh_varnan] = new_var[wh_varnan]
    variance[wh_numer0] = new_var[wh_numer0]

    return variance


def calc_num_seg(gdq, imshape, n_int, nreads):
    """
    Extended Summary
    ----------------
    Calculate the maximum number of segments that will be fit within an
    integration, calculated over all pixels and all integrations. This
    value will be used to allocate arrays used for the optional output
    product, and is based on the locations of saturated pixels and cosmic
    ray-affected pixels in all of the ramps.

    Parameters
    ----------
    gdq: float, 3D array
        cube of GROUPDQ array for a data

    imshape: (int, int) tuple
        shape of 2D image

    n_int: int
        total number of integrations in data set

    nreads: int
        number of reads in an integration

    Return:
    -------
    num_seg.max(): int
        maximum number of segments
    """

    # Initialize pixel array for number of segments for each integration
    num_seg = np.zeros((n_int,) + imshape, dtype=np.int8)

    for nint in range(n_int):  # loop over integrations

        gdq_int = gdq[nint, :, :, :] # this integration's GROUPDQ cube

        # each pixel's number of segments for current integration:
        num_seg_int = num_seg[nint, :, :]

        # Initialize pixel array to designate where previous read is good
        #   (which means 'not a CR'); here the (fictitious) read previous
        #   to the 1st read is set to False, because the number of segments
        #   will be incremented when a second consecutive True previous read
        #   after a False is encounterd.
        prev_good = np.zeros(imshape, dtype=np.bool)
        prev_good[:] = False

        # Initialize pixel array designating which pixels to still process;
        #   1 means still to do
        pix_do = np.ones(imshape, dtype=np.int)

        # Find pixels having all good reads
        wh_good = np.where(gdq_int.sum(axis=0) == 0)
        num_seg_int[wh_good] = 1 # 1 segment containing all reads
        pix_do[wh_good] = 0 # these all good pix are done
        prev_good[wh_good] = True

        # Find pixels with SAT on 1st read; these will have 0 segments
        #   fit in this integration,
        wh_e_sat = np.where(np.bitwise_and(gdq_int[0, :, :],
                                             dqflags.group['SATURATED']) != 0)

        pix_do[wh_e_sat] = 0  # flag these early sat pix as done

        # Assign vectors to for loop that follows
        num_seg_1d = num_seg_int.ravel()
        pix_do_1d = pix_do.ravel()
        prev_1d = prev_good.ravel()

        # Loop over all reads, selecting pixels not yet done for further
        #   processing
        for rd in range(nreads):
            sl_rd = gdq_int[rd, :, :].ravel()

            # Find saturated pixels and flag as done
            wh_sat = np.where((pix_do_1d == 1) &
                    (np.bitwise_and(sl_rd, dqflags.group['SATURATED']) != 0))

            pix_do_1d[wh_sat] = 0

            # For cr-affected pixels following good read, flag 'previous' array
            #   as False so next good read will signify a segment end
            wh_cr = np.where((pix_do_1d == 1) &
                 (np.bitwise_and(sl_rd, dqflags.group['JUMP_DET']) != 0) &
                 (prev_1d == True))
            prev_1d[wh_cr] = False

            # For good reads immediately following a cosmic ray, increment the
            #   number of segments and flag the 'previous' read as 'True'. This
            #   doesn't apply if this is the 1st read, as the segment may not be
            #   long enough
            wh_ok = np.where((pix_do_1d == 1) & (sl_rd == 0) &
                             (prev_1d == False) & (rd > 0))

            num_seg_1d[wh_ok] += 1 # this also increments num_seg_int
            prev_1d[wh_ok] = True

        num_seg[nint, :, :] = num_seg_int.reshape(imshape)

    # log and return global maximum
    log.debug('Maximum number of segments per integration:%d' % (num_seg.max()))

    return num_seg.max()


def calc_unwtd_sums(data_masked, xvalues):
    """
    Short Summary
    -------------
    Calculate the sums needed to determine the slope and intercept (and sigma
    of each) using an unweighted fit.

    Parameters
    ----------
    data_masked: float, 2D array
        masked values for all pixels in data section

    xvalues: int, 1D array
        indices of valid pixel values for all reads

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


def calc_opt_sums(rn_sect, gain_sect, data_masked, mask_2d, xvalues):
    """
    Short Summary
    -------------
    Calculate the sums needed to determine the slope and intercept (and
    sigma of each) using the optimal weights.  For each good pixel's
    semiramp, from the initial and final indices and the corresponding
    number of counts, calculate the SNR. From the SNR, calculate the
    weighting exponent using the formulation by Fixsen (Fixsen et al, PASP,
    112, 1350). Using this exponent and the gain and the readnoise, the weights
    are calculated from which the sums are calculated.

    Parameters
    ----------
    rn_sect: float, 2D array
        read noise values for all pixels in data section

    gain_sect: float, 2D array
        gain values for all pixels in data section

    data_masked: float, 2D array
        masked values for all pixels in data section

    mask_2d: int, 2D array
        delineates which channels to fit for each pixel

    xvalues: int, 2D array
        indices of valid pixel values for all reads

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
        sum of product of data and weight

    xvalues: int, 2D array
        rolled up indices of valid pixel values for all reads
    """

    # Return 'empty' sums if there is no more data to fit
    if (data_masked.size == 0):
        return np.array([]), np.array([]), np.array([]), np.array([]),\
               np.array([]), np.array([])

    # get initial and final valid reads for each pixel
    fnz = np.argmax(mask_2d, axis=0) # lowest read index that is True

    # For those pixels that are all False, set to sentinel value of -1:
    fnz[mask_2d.sum(axis=0) == 0] = -1

    mask_2d_sum = mask_2d.sum(axis=0)   # number of valid reads/pixel
    ind_lastnz = fnz + mask_2d_sum - 1
    data_zero = data_masked[fnz, range(data_masked.shape[1])]
    data_final = data_masked[(ind_lastnz), range(data_masked.shape[1])]
    data_diff = data_final - data_zero

    ind_lastnz = 0
    data_zero = 0

   # Use the readnoise and gain for good pixels only
    rn_2_r = (rn_sect * rn_sect).ravel()
    rn_2_r = rn_2_r[fnz]

    gain_sect_r = gain_sect.ravel()
    gain_sect = 0
    gain_sect_r = gain_sect_r[fnz]

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

    gain_sect_r = 0
    numer_ir = 0
    data_diff = 0
    sigma_ir = 0

    power_wt_r = calc_power(snr)  # get the weighting exponent for this SNR

   # Make array of number of good reads, and exponents for each pixel
    num_nz = (data_masked != 0.).sum(0) # number of nonzero reads per pixel
    nrd_data_a = num_nz.copy()
    num_nz = 0

    nrd_prime = (nrd_data_a - 1) / 2.
    nrd_data_a = 0

    # Calculate inverse read noise^2 for use in weights
    invrdns2_r = (1. / (rn_sect * rn_sect)).ravel()
    rn_sect = 0
    invrdns2_r = invrdns2_r[fnz]
    fnz = 0

    # Set optimal weights for each read of each pixel;
    #    for all pixels at once, loop over the read
    wt_h = np.zeros(data_masked.shape, dtype=np.float32)

    for jj_rd in range(data_masked.shape[0]):
        wt_h[jj_rd, :] = \
              abs((abs(jj_rd - nrd_prime) / nrd_prime) ** power_wt_r) * invrdns2_r

    wt_h[np.isnan(wt_h)] = 0.
    wt_h[np.isinf(wt_h)] = 0.

    nrd_prime = 0
    power_wt_r = 0
    invrdns2_r = 0

    wh_m2d_f = (mask_2d[0, :] == False) # where initial read is False

    # For all pixels, 'roll' up the leading zeros such that the 0th read of
    #  each pixel is the lowest nonzero read for that pixel
    while (wh_m2d_f.sum() > 0):
        data_masked[:, wh_m2d_f] = np.roll(data_masked[:, wh_m2d_f], -1, axis=0)
        mask_2d[:, wh_m2d_f] = np.roll(mask_2d[:, wh_m2d_f], -1, axis=0)
        xvalues[:, wh_m2d_f] = np.roll(xvalues[:, wh_m2d_f], -1, axis=0)
        wh_m2d_f = (mask_2d[0, :] == False)

    wt_h[data_masked == 0.] = 0. # 0 the weight for pixels with all reads=0

    # Create sums
    nreads_wtd = wt_h.sum(axis=0)  # sum of each pixel's weights
    sumx = (xvalues * wt_h).sum(axis=0)
    sumxx = (xvalues**2 * wt_h).sum(axis=0)
    sumy = (np.reshape((data_masked * wt_h).sum(axis=0), sumx.shape))
    sumxy = (xvalues * wt_h * np.reshape(data_masked, xvalues.shape)).sum(axis=0)

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

    log.debug('The number of pixels having insufficient data'),
    log.debug('due to excessive CRs or saturation %d:' % (len(wh_c_0[0])))
    log.debug('Count rates - min, mean, max, std: %f, %f, %f, %f'
               % (c_rates.min(), c_rates.mean(), c_rates.max(), c_rates.std()))
