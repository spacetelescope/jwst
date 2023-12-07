#! /usr/bin/env python

import numpy as np

from stcal.ramp_fitting import ramp_fit
from stcal.ramp_fitting import utils

from stcal.ramp_fitting.utils import LARGE_VARIANCE
from stcal.ramp_fitting.utils import LARGE_VARIANCE_THRESHOLD

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from ..stpipe import Step

from ..lib import reffile_utils

import logging
import copy
import warnings
import multiprocessing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

multiprocessing.set_start_method('forkserver', force=True)


__all__ = ["RampFitStep"]


def get_reference_file_subarrays(model, readnoise_model, gain_model, nframes):
    """
    Get readnoise array for calculation of variance of noiseless ramps, and
    the gain array in case optimal weighting is to be done. The returned
    readnoise has been multiplied by the gain.

    Parameters
    ----------
    model : data model
        input data model, assumed to be of type RampModel

    readnoise_model : instance of data Model
        readnoise for all pixels

    gain_model : instance of gain Model
        gain for all pixels

    nframes : int
        number of frames averaged per group; from the NFRAMES keyword. Does
        not contain the groupgap.

    Returns
    -------
    readnoise_2d : float, 2D array
        readnoise subarray

    gain_2d : float, 2D array
        gain subarray
    """
    if reffile_utils.ref_matches_sci(model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info('Extracting gain subarray to match science data')
        gain_2d = reffile_utils.get_subarray_data(model, gain_model)

    if reffile_utils.ref_matches_sci(model, readnoise_model):
        readnoise_2d = readnoise_model.data.copy()
    else:
        log.info('Extracting readnoise subarray to match science data')
        readnoise_2d = reffile_utils.get_subarray_data(model, readnoise_model)

    return readnoise_2d, gain_2d


def create_image_model(input_model, image_info):
    """
    Creates an ImageModel from the computed arrays from ramp_fit.

    Parameter
    ---------
    input_model: RampModel
        Input RampModel for which the output ImageModel is created.

    image_info: tuple
        The ramp fitting arrays needed for the ImageModel.

    Return
    ---------
    out_model: ImageModel
        The output ImageModel to be returned from the ramp fit step.
    """
    data, dq, var_poisson, var_rnoise, err = image_info

    # Create output datamodel
    out_model = datamodels.ImageModel(data.shape)

    # ... and add all keys from input
    out_model.update(input_model)

    # Populate with output arrays
    out_model.data = data
    out_model.dq = dq
    out_model.var_poisson = var_poisson
    out_model.var_rnoise = var_rnoise
    out_model.err = err

    return out_model


def create_integration_model(input_model, integ_info, int_times):
    """
    Creates an ImageModel from the computed arrays from ramp_fit.

    Parameter
    ---------
    input_model : RampModel
        Input RampModel for which the output CubeModel is created.

    integ_info: tuple
        The ramp fitting arrays needed for the CubeModel for each integration.

    int_times : astropy.io.fits.fitsrec.FITS_rec or None
        Integration times.

    Return
    ---------
    int_model : CubeModel
        The output CubeModel to be returned from the ramp fit step.
    """
    data, dq, var_poisson, var_rnoise, err = integ_info
    int_model = datamodels.CubeModel(
        data=np.zeros(data.shape, dtype=np.float32),
        dq=np.zeros(data.shape, dtype=np.uint32),
        var_poisson=np.zeros(data.shape, dtype=np.float32),
        var_rnoise=np.zeros(data.shape, dtype=np.float32),
        err=np.zeros(data.shape, dtype=np.float32))
    int_model.update(input_model)  # ... and add all keys from input

    int_model.data = data
    int_model.dq = dq
    int_model.var_poisson = var_poisson
    int_model.var_rnoise = var_rnoise
    int_model.err = err
    int_model.int_times = int_times

    return int_model


def create_optional_results_model(input_model, opt_info):
    """
    Creates an ImageModel from the computed arrays from ramp_fit.

    Parameter
    ---------
    input_model: ~jwst.datamodels.RampModel

    opt_info: tuple
        The ramp fitting arrays needed for the RampFitOutputModel.

    Return
    ---------
    opt_model: RampFitOutputModel
        The optional RampFitOutputModel to be returned from the ramp fit step.
    """
    (slope, sigslope, var_poisson, var_rnoise,
        yint, sigyint, pedestal, weights, crmag) = opt_info
    opt_model = datamodels.RampFitOutputModel(
        slope=slope,
        sigslope=sigslope,
        var_poisson=var_poisson,
        var_rnoise=var_rnoise,
        yint=yint,
        sigyint=sigyint,
        pedestal=pedestal,
        weights=weights,
        crmag=crmag)

    opt_model.meta.filename = input_model.meta.filename
    opt_model.update(input_model)  # ... and add all keys from input

    return opt_model


def compute_RN_variances(groupdq, readnoise_2d, gain_2d, group_time):
    """
    Compute the variances due to the readnoise for all integrations.

    Parameters
    ----------
    groupdq : ndarray
        The group data quality array for the exposure, 4-D flag.
        For groups that have been flagged as both CHARGELOSS and
        DO_NOT_USE, both flags have been reset.

    readnoise_2d : ndarray
        readnoise values for all pixels in the image, 2-D float

    gain_2d : ndarray
        gain values for all pixels in the image, 2-D float

    group_time : float
        Time increment between groups, in seconds.

    Returns
    -------
    var_r2 : ndarray
        Image of integration-specific values for the slope variance due to
        readnoise only, 2-D float

    var_r3 : ndarray
        Cube of integration-specific values for the slope variance due to
        readnoise only, 3-D float

    var_r4 : ndarray
        Hypercube of segment- and integration-specific values for the slope
        variance due to read noise only, 4-D float
    """
    nint, ngroups, nrows, ncols = groupdq.shape

    imshape = (nrows, ncols)
    cubeshape = (ngroups,) + imshape

    segs_4 = np.zeros((nint,) + (ngroups,) + imshape, dtype=np.uint8)
    var_r4 = np.zeros((nint,) + (ngroups,) + imshape, dtype=np.float32) + LARGE_VARIANCE
    var_r3 = np.zeros((nint,) + imshape, dtype=np.float32) + LARGE_VARIANCE
    s_inv_var_r3 = np.zeros((nint,) + imshape, dtype=np.float32)

    max_seg = 0  # initialize maximum number of segments in a ramp

    # Loop over data integrations
    for num_int in range(nint):
        # Loop over data sections
        for rlo in range(0, cubeshape[1], nrows):
            rhi = rlo + nrows

            if rhi > cubeshape[1]:
                rhi = cubeshape[1]

            gdq_sect = groupdq[num_int, :, rlo:rhi, :]
            rn_sect = readnoise_2d[rlo:rhi, :]
            gain_sect = gain_2d[rlo:rhi, :]

            # For each data section, calculate values of segment lengths and
            #   quantities to calculate variances.
            den_r3, num_r3, segs_beg_3, max_seg_i = calc_segs(rn_sect, gdq_sect, group_time)
            max_seg = max(max_seg, max_seg_i)
            segs_4[num_int, :, rlo:rhi, :] = segs_beg_3

            # Find the segment variance due to read noise and convert back to DN
            var_r4[num_int, :, rlo:rhi, :] = num_r3 * den_r3 / gain_sect**2

            del den_r3, num_r3, segs_beg_3
            del gdq_sect
            del rn_sect
            del gain_sect

        # Zero out entries for non-existing segments in this integration
        var_r4[num_int, :, :, :] *= (segs_4[num_int, :, :, :] > 0)

        # Correct for non-positive values to ensure negligible conntribution
        var_r4[var_r4 <= 0.] = LARGE_VARIANCE

        # The sums of inverses of the variances are needed for later
        #   variance calculations.
        s_inv_var_r3[num_int, :, :] = (1. / var_r4[num_int, :, :, :]).sum(axis=0)
        var_r3[num_int, :, :] = 1. / s_inv_var_r3[num_int, :, :]

    var_r4 *= (segs_4[:, :, :, :] > 0)

    # Truncate var_r4 to include only existing segments
    var_r4 = var_r4[:, :max_seg, :, :]

    var_r3[var_r3 > LARGE_VARIANCE_THRESHOLD] = 0.  # Zero out large variances due to reset
    var_r2 = 1 / (s_inv_var_r3.sum(axis=0))
    var_r2[var_r2 > LARGE_VARIANCE_THRESHOLD] = 0.
    var_r2[np.isnan(var_r2)] = 0.

    return var_r2, var_r3, var_r4


def calc_segs(rn_sect, gdq_sect, group_time):
    """
    Calculate several quantities needed for the readnoise variance, in the
    data section.

    Parameters
    ----------
    rn_sect : ndarray
        readnoise values for all pixels in the data section, 2-D float

    gdq_sect : ndarray
        gain values for all pixels in the data section, 2-D float

    group_time : float
        Time increment between groups, in seconds.

    Returns
    -------
    den_r3 : ndarray
        for a given integration, the reciprocal of the denominator of the
        segment-specific variance of the segment's slope due to read noise, 3-D float

    num_r3 : ndarray
        numerator of the segment-specific variance of the segment's slope
        due to read noise, 3-D float

    segs_beg_3 : ndarray
        lengths of segments for all pixels in the given data section and
        integration, 3-D

    max_seg : int
        maximum number of segments in a ramp
    """
    (ngroups, asize2, asize1) = gdq_sect.shape
    npix = asize1 * asize2
    imshape = (asize2, asize1)
    gdq_2d = gdq_sect[:, :, :].reshape((ngroups, npix))
    segs = np.zeros((ngroups, npix), dtype=np.int32)
    sr_index = np.zeros(npix, dtype=np.uint8)

    i_read = 0
    while i_read < ngroups:
        gdq_1d = gdq_2d[i_read, :]
        wh_good = np.where(gdq_1d == dqflags.group['GOOD'])

        # For good groups, increment those pixels' segments' lengths
        if len(wh_good[0]) > 0:
            segs[sr_index[wh_good], wh_good] += 1
        del wh_good

        # Locate any CRs ...
        wh_cr = np.where(gdq_1d.astype(np.int32) & dqflags.group['JUMP_DET'] > 0)
        del gdq_1d

        # ... (but not on final read): increment the segment number
        if len(wh_cr[0]) > 0 and (i_read < ngroups - 1):
            sr_index[wh_cr[0]] += 1
            segs[sr_index[wh_cr], wh_cr] += 1
        del wh_cr

        i_read += 1

    segs = segs.astype(np.uint8)
    segs_beg_3 = segs.reshape(ngroups, imshape[0], imshape[1])
    segs_beg_3 = utils.remove_bad_singles(segs_beg_3)

    # For a segment, the variance due to readnoise noise
    # = 12 * readnoise**2 /(nreads_seg**3. - nreads_seg)/(tgroup **2.)
    num_r3 = 12. * (rn_sect / group_time)**2.  # always >0

    # Reshape for every group, every pixel in section
    num_r3 = np.dstack([num_r3] * ngroups)
    num_r3 = np.transpose(num_r3, (2, 0, 1))

    # Denominator den_r3 = 1./(segs_beg_3 **3.-segs_beg_3). The minimum number
    #   of allowed groups is 2, which will apply if there is actually only 1
    #   group; in this case den_r3 = 1/6. This covers the case in which there is
    #   only one good group at the beginning of the integration, so it will be
    #   be compared to the plane of (near) zeros resulting from the reset. For
    #   longer segments, this value is overwritten below.
    den_r3 = num_r3.copy() * 0. + 1. / 6

    wh_seg_pos = np.where(segs_beg_3 > 1)

    # Suppress, then, re-enable harmless arithmetic warnings, as NaN will be
    #   checked for and handled later
    warnings.filterwarnings("ignore", ".*invalid value.*", RuntimeWarning)
    warnings.filterwarnings("ignore", ".*divide by zero.*", RuntimeWarning)
    # overwrite where segs>1
    den_r3[wh_seg_pos] = 1. / (segs_beg_3[wh_seg_pos] ** 3. - segs_beg_3[wh_seg_pos])
    warnings.resetwarnings()

    # calculate max_seg for this integ and data section
    max_seg = (np.count_nonzero(segs_beg_3, axis=0)).max()

    return den_r3, num_r3, segs_beg_3, max_seg


class RampFitStep(Step):

    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    class_alias = "ramp_fit"

    spec = """
        int_name = string(default='')
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        suppress_one_group = boolean(default=True)  # Suppress saturated ramps with good 0th group
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
    """

    # Prior to 04/26/17, the following were also in the spec above:
    #      algorithm = option('OLS', 'GLS', default='OLS') # 'OLS' or 'GLS'
    #      weighting = option('unweighted', 'optimal', default='unweighted') \
    #      # 'unweighted' or 'optimal'
    # As of 04/26/17, the only allowed algorithm is 'ols', and the
    #      only allowed weighting is 'optimal'.

    algorithm = 'ols'      # Only algorithm allowed for Build 7.1
#    algorithm = 'gls'       # 032520

    weighting = 'optimal'  # Only weighting allowed for Build 7.1

    reference_file_types = ['readnoise', 'gain']

    def process(self, input):

        with datamodels.RampModel(input) as input_model:
            max_cores = self.maximum_cores
            readnoise_filename = self.get_reference_file(input_model, 'readnoise')
            gain_filename = self.get_reference_file(input_model, 'gain')

            log.info('Using READNOISE reference file: %s', readnoise_filename)
            log.info('Using GAIN reference file: %s', gain_filename)

            with datamodels.ReadnoiseModel(readnoise_filename) as readnoise_model, \
                    datamodels.GainModel(gain_filename) as gain_model:

                # Try to retrieve the gain factor from the gain reference file.
                # If found, store it in the science model meta data, so that it's
                # available later in the gain_scale step, which avoids having to
                # load the gain ref file again in that step.
                if gain_model.meta.exposure.gain_factor is not None:
                    input_model.meta.exposure.gain_factor = gain_model.meta.exposure.gain_factor

                # Get gain arrays, subarrays if desired.
                frames_per_group = input_model.meta.exposure.nframes
                readnoise_2d, gain_2d = get_reference_file_subarrays(
                    input_model, readnoise_model, gain_model, frames_per_group)

            log.info('Using algorithm = %s' % self.algorithm)
            log.info('Using weighting = %s' % self.weighting)

            buffsize = ramp_fit.BUFSIZE
            if self.algorithm == "GLS":
                buffsize //= 10

            int_times = input_model.int_times

            # Before the ramp_fit() call, copy the input model ("_W" for weighting)
            # for later reconstruction of the fitting array tuples.
            input_model_W = copy.copy(input_model)

            # Run ramp_fit(), ignoring all DO_NOT_USE groups, and return the
            # ramp fitting arrays for the ImageModel, the CubeModel, and the
            # RampFitOutputModel.
            image_info, integ_info, opt_info, gls_opt_model = ramp_fit.ramp_fit(
                input_model, buffsize, self.save_opt, readnoise_2d, gain_2d,
                self.algorithm, self.weighting, max_cores, dqflags.pixel,
                suppress_one_group=self.suppress_one_group)

            # Create a gdq to modify if there are charge_migrated groups
            gdq = input_model_W.groupdq.copy()

            # Locate groups where that are flagged with CHARGELOSS
            wh_chargeloss = np.where(np.bitwise_and(gdq.astype(np.uint32), dqflags.group['CHARGELOSS']))

            if len(wh_chargeloss[0]) > 0:
                # Unflag groups flagged as both CHARGELOSS and DO_NOT_USE
                gdq[wh_chargeloss] -= (dqflags.group['DO_NOT_USE'] + dqflags.group['CHARGELOSS'])

                # Flag SATURATED groups as DO_NOT_USE for later segment determination
                where_sat = np.where(np.bitwise_and(gdq, dqflags.group['SATURATED']))
                gdq[where_sat] = np.bitwise_or(gdq[where_sat], dqflags.group['DO_NOT_USE'])

                # Get group_time for readnoise variance calculation
                group_time = input_model.meta.exposure.group_time

                # Using the modified GROUPDQ array, create new readnoise variance arrays
                image_var_RN, integ_var_RN, opt_var_RN = \
                    compute_RN_variances(gdq, readnoise_2d, gain_2d, group_time)

                # Create new ramp fitting array tuples, by inserting the new
                # readnoise variances into copies of the original ramp fitting
                # tuples.
                image_info_new, integ_info_new = None, None
                if image_info is not None and image_var_RN is not None:
                    image_info_new = (image_info[0], image_info[1], image_info[2], image_var_RN, image_info[4])

                if integ_info is not None and integ_var_RN is not None:
                    integ_info_new = (integ_info[0], integ_info[1], integ_info[2], integ_var_RN, integ_info[4])

                image_info = image_info_new
                integ_info = integ_info_new

                opt_info_new = None
                if opt_info is not None and opt_var_RN is not None:
                    opt_info_new = (opt_info[0], opt_info[1], opt_info[2], opt_var_RN,
                                    opt_info[4], opt_info[5], opt_info[6], opt_info[7], opt_info[8])

                opt_info = opt_info_new

        # Save the OLS optional fit product, if it exists.
        if opt_info is not None:
            opt_model = create_optional_results_model(input_model, opt_info)
            self.save_model(opt_model, 'fitopt', output_file=self.opt_name)
        '''
        # GLS removed from code, since it's not implemented right now.
        # Save the GLS optional fit product, if it exists
        if gls_opt_model is not None:
            self.save_model(
                gls_opt_model, 'fitoptgls', output_file=self.opt_name
            )
        '''

        out_model, int_model = None, None
        # Create models from possibly updated info
        if image_info is not None and integ_info is not None:
            out_model = create_image_model(input_model, image_info)
            out_model.meta.bunit_data = 'DN/s'
            out_model.meta.bunit_err = 'DN/s'
            out_model.meta.cal_step.ramp_fit = 'COMPLETE'
            if ((input_model.meta.exposure.type in ['NRS_IFU', 'MIR_MRS']) or
                    (input_model.meta.exposure.type in ['NRS_AUTOWAVE', 'NRS_LAMP'] and
                     input_model.meta.instrument.lamp_mode == 'IFU')):

                out_model = datamodels.IFUImageModel(out_model)

            int_model = create_integration_model(input_model, integ_info, int_times)
            int_model.meta.bunit_data = 'DN/s'
            int_model.meta.bunit_err = 'DN/s'
            int_model.meta.cal_step.ramp_fit = 'COMPLETE'

        return out_model, int_model
