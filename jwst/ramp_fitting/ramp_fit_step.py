#! /usr/bin/env python
import numpy as np

from stcal.ramp_fitting import ramp_fit
from stcal.ramp_fitting import utils

from stcal.ramp_fitting.likely_fit import LIKELY_MIN_NGROUPS

from stcal.ramp_fitting.utils import LARGE_VARIANCE
from stcal.ramp_fitting.utils import LARGE_VARIANCE_THRESHOLD

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags


from ..stpipe import Step

from ..lib import reffile_utils

import logging
import warnings

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["RampFitStep"]


def get_reference_file_subarrays(model, readnoise_model, gain_model):
    """
    Get read noise and gain reference arrays.

    Get readnoise array for calculation of variance of noiseless ramps, and
    the gain array in case optimal weighting is to be done. The returned
    readnoise has been multiplied by the gain.

    Parameters
    ----------
    model : data model
        Input data model, assumed to be of type RampModel.
    readnoise_model : instance of data Model
        Readnoise for all pixels.
    gain_model : instance of gain Model
        Gain for all pixels.

    Returns
    -------
    readnoise_2d : float, 2D array
        Readnoise subarray
    gain_2d : float, 2D array
        Gain subarray
    """
    if reffile_utils.ref_matches_sci(model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info("Extracting gain subarray to match science data")
        gain_2d = reffile_utils.get_subarray_model(model, gain_model).data

    if reffile_utils.ref_matches_sci(model, readnoise_model):
        readnoise_2d = readnoise_model.data
    else:
        log.info("Extracting readnoise subarray to match science data")
        readnoise_2d = reffile_utils.get_subarray_model(model, readnoise_model).data

    return readnoise_2d, gain_2d


def create_image_model(input_model, image_info):
    """
    Create an ImageModel from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : RampModel
        Input RampModel for which the output ImageModel is created.
    image_info : tuple
        The ramp fitting arrays needed for the ImageModel.

    Returns
    -------
    out_model : ImageModel
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
    Create an ImageModel from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : RampModel
        Input RampModel for which the output CubeModel is created.
    integ_info : tuple
        The ramp fitting arrays needed for the CubeModel for each integration.
    int_times : astropy.io.fits.fitsrec.FITS_rec or None
        Integration times.

    Returns
    -------
    int_model : CubeModel
        The output CubeModel to be returned from the ramp fit step.
    """
    data, dq, var_poisson, var_rnoise, err = integ_info
    int_model = datamodels.CubeModel(
        data=np.zeros(data.shape, dtype=np.float32),
        dq=np.zeros(data.shape, dtype=np.uint32),
        var_poisson=np.zeros(data.shape, dtype=np.float32),
        var_rnoise=np.zeros(data.shape, dtype=np.float32),
        err=np.zeros(data.shape, dtype=np.float32),
    )
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
    Create an ImageModel from the computed arrays from ramp_fit.

    Parameters
    ----------
    input_model : RampModel
        The input ramp model used to create the optional results product.
    opt_info : tuple
        The ramp fitting arrays needed for the RampFitOutputModel.

    Returns
    -------
    opt_model : RampFitOutputModel
        The optional RampFitOutputModel to be returned from the ramp fit step.
    """
    opt_model = datamodels.RampFitOutputModel(
        slope=opt_info[0],
        sigslope=opt_info[1],
        var_poisson=opt_info[2],
        var_rnoise=opt_info[3],
        yint=opt_info[4],
        sigyint=opt_info[5],
        pedestal=opt_info[6],
        weights=opt_info[7],
        crmag=opt_info[8],
    )

    opt_model.meta.filename = input_model.meta.filename
    opt_model.update(input_model)  # ... and add all keys from input

    return opt_model


def compute_rn_variances(groupdq, readnoise_2d, gain_2d, group_time):
    """
    Compute the variances due to the readnoise for all integrations.

    Parameters
    ----------
    groupdq : ndarray
        The group data quality array for the exposure, 4-D flag.
        For groups that have been flagged as both CHARGELOSS and
        DO_NOT_USE, both flags have been reset.
    readnoise_2d : ndarray
        Readnoise values for all pixels in the image, 2-D float.
    gain_2d : ndarray
        Gain values for all pixels in the image, 2-D float.
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

    segs_4 = np.zeros((nint,) + (ngroups,) + imshape, dtype=np.uint16)
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
        var_r4[num_int, :, :, :] *= segs_4[num_int, :, :, :] > 0

        # Correct for non-positive values to ensure negligible conntribution
        var_r4[var_r4 <= 0.0] = LARGE_VARIANCE

        # The sums of inverses of the variances are needed for later
        #   variance calculations.
        s_inv_var_r3[num_int, :, :] = (1.0 / var_r4[num_int, :, :, :]).sum(axis=0)
        var_r3[num_int, :, :] = 1.0 / s_inv_var_r3[num_int, :, :]

    var_r4 *= segs_4[:, :, :, :] > 0

    # Truncate var_r4 to include only existing segments
    var_r4 = var_r4[:, :max_seg, :, :]

    # Zero out large variances due to reset
    var_r3[var_r3 > LARGE_VARIANCE_THRESHOLD] = 0.0
    var_r2 = 1 / (s_inv_var_r3.sum(axis=0))
    var_r2[var_r2 > LARGE_VARIANCE_THRESHOLD] = 0.0
    var_r2[np.isnan(var_r2)] = 0.0

    return var_r2, var_r3, var_r4


def calc_segs(rn_sect, gdq_sect, group_time):
    """
    Calculate quantities needed for the readnoise variance.

    Parameters
    ----------
    rn_sect : ndarray
        Readnoise values for all pixels in the data section, 2-D float.
    gdq_sect : ndarray
        Gain values for all pixels in the data section, 2-D float.
    group_time : float
        Time increment between groups, in seconds.

    Returns
    -------
    den_r3 : ndarray
        For a given integration, the reciprocal of the denominator
        of the segment-specific variance of the segment's slope due
        to read noise, 3-D float.
    num_r3 : ndarray
        Numerator of the segment-specific variance of the segment's slope
        due to read noise, 3-D float.
    segs_beg_3 : ndarray
        Lengths of segments for all pixels in the given data section and
        integration, 3-D.
    max_seg : int
        Maximum number of segments in a ramp.
    """
    (ngroups, asize2, asize1) = gdq_sect.shape
    npix = asize1 * asize2
    imshape = (asize2, asize1)
    gdq_2d = gdq_sect[:, :, :].reshape((ngroups, npix))
    segs = np.zeros((ngroups, npix), dtype=np.int32)
    sr_index = np.zeros(npix, dtype=np.uint16)

    i_read = 0
    while i_read < ngroups:
        gdq_1d = gdq_2d[i_read, :]
        wh_good = np.where(gdq_1d == dqflags.group["GOOD"])

        # For good groups, increment those pixels' segments' lengths
        if len(wh_good[0]) > 0:
            segs[sr_index[wh_good], wh_good] += 1
        del wh_good

        # Locate any CRs ...
        wh_cr = np.where(gdq_1d.astype(np.int32) & dqflags.group["JUMP_DET"] > 0)
        del gdq_1d

        # ... (but not on final read): increment the segment number
        if len(wh_cr[0]) > 0 and (i_read < ngroups - 1):
            sr_index[wh_cr[0]] += 1
            segs[sr_index[wh_cr], wh_cr] += 1
        del wh_cr

        i_read += 1

    segs = segs.astype(np.uint16)
    segs_beg_3 = segs.reshape(ngroups, imshape[0], imshape[1])
    segs_beg_3 = utils.remove_bad_singles(segs_beg_3)

    # For a segment, the variance due to readnoise noise
    # = 12 * readnoise**2 /(nreads_seg**3. - nreads_seg)/(tgroup **2.)
    num_r3 = 12.0 * (rn_sect / group_time) ** 2.0  # always >0

    # Reshape for every group, every pixel in section
    num_r3 = np.dstack([num_r3] * ngroups)
    num_r3 = np.transpose(num_r3, (2, 0, 1))

    # Denominator den_r3 = 1./(segs_beg_3 **3.-segs_beg_3). The minimum number
    #   of allowed groups is 2, which will apply if there is actually only 1
    #   group; in this case den_r3 = 1/6. This covers the case in which there is
    #   only one good group at the beginning of the integration, so it will be
    #   be compared to the plane of (near) zeros resulting from the reset. For
    #   longer segments, this value is overwritten below.
    den_r3 = num_r3.copy() * 0.0 + 1.0 / 6

    wh_seg_pos = np.where(segs_beg_3 > 1)

    # Suppress, then, re-enable harmless arithmetic warnings, as NaN will be
    #   checked for and handled later
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", ".*invalid value.*", RuntimeWarning)
        warnings.filterwarnings("ignore", ".*divide by zero.*", RuntimeWarning)
        # overwrite where segs>1
        den = segs_beg_3[wh_seg_pos] ** 3.0 - segs_beg_3[wh_seg_pos]
        den_r3[wh_seg_pos] = 1.0 / den

    # calculate max_seg for this integ and data section
    max_seg = (np.count_nonzero(segs_beg_3, axis=0)).max()

    return den_r3, num_r3, segs_beg_3, max_seg


class RampFitStep(Step):
    """Fit line to determine the value of mean rate counts vs. time."""

    class_alias = "ramp_fit"

    spec = """
        algorithm = option('OLS', 'OLS_C', 'LIKELY', default='OLS_C') # 'OLS' and 'OLS_C' use the same underlying algorithm, but OLS_C is implemented in C
        int_name = string(default='')
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        suppress_one_group = boolean(default=True)  # Suppress saturated ramps with good 0th group
        firstgroup = integer(default=None)   # Ignore groups before this one (zero indexed)
        lastgroup = integer(default=None)   # Ignore groups after this one (zero indexed)
        maximum_cores = string(default='1') # cores for multiprocessing. Can be an integer, 'half', 'quarter', or 'all'
    """  # noqa: E501

    # Prior to 04/26/17, the following were also in the spec above:
    #      algorithm = option('OLS', 'GLS', default='OLS') # 'OLS' or 'GLS'
    #      weighting = option('unweighted', 'optimal', default='unweighted') \
    #      # 'unweighted' or 'optimal'
    # As of 04/26/17, the only allowed algorithm is 'ols', and the
    #      only allowed weighting is 'optimal'.

    # algorithm = 'ols'      # Only algorithm allowed for Build 7.1
    # algorithm = 'gls'       # 032520

    weighting = "optimal"  # Only weighting allowed for Build 7.1

    reference_file_types = ["readnoise", "gain"]

    def process(self, step_input):
        """
        Fit ramps using the specified ramp fitting algorithm.

        Parameters
        ----------
        step_input : RampModel
            The input ramp model to fit the ramps.

        Returns
        -------
        out_model : ImageModel
            The output 2-D image model with the fit ramps.

        int_model : CubeModel
            The output 3-D image model with the fit ramps for each integration.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # Work on a copy
            result = input_model.copy()

            max_cores = self.maximum_cores
            readnoise_filename = self.get_reference_file(result, "readnoise")
            gain_filename = self.get_reference_file(result, "gain")

            ngroups = input_model.data.shape[1]
            if self.algorithm.upper() == "LIKELY" and ngroups < LIKELY_MIN_NGROUPS:
                log.info(
                    f"When selecting the LIKELY ramp fitting algorithm the"
                    f" ngroups needs to be a minimum of {LIKELY_MIN_NGROUPS},"
                    f" but ngroups = {ngroups}.  Due to this, the ramp fitting algorithm"
                    f" is being changed to OLS_C"
                )
                self.algorithm = "OLS_C"

            log.info(f"Using READNOISE reference file: {readnoise_filename}")
            log.info(f"Using GAIN reference file: {gain_filename}")

            with (
                datamodels.ReadnoiseModel(readnoise_filename) as readnoise_model,
                datamodels.GainModel(gain_filename) as gain_model,
            ):
                # Try to retrieve the gain factor from the gain reference file.
                # If found, store it in the science model meta data, so that it's
                # available later in the gain_scale step, which avoids having to
                # load the gain ref file again in that step.
                if gain_model.meta.exposure.gain_factor is not None:
                    result.meta.exposure.gain_factor = gain_model.meta.exposure.gain_factor

                # Get gain arrays, subarrays if desired.
                readnoise_2d, gain_2d = get_reference_file_subarrays(
                    result, readnoise_model, gain_model
                )

            log.info(f"Using algorithm = {self.algorithm}")
            log.info(f"Using weighting = {self.weighting}")

            buffsize = ramp_fit.BUFSIZE
            if self.algorithm == "GLS":
                buffsize //= 10

            int_times = result.int_times

            # Set the DO_NOT_USE bit in the groupdq values for groups before firstgroup
            # and groups after lastgroup
            firstgroup = self.firstgroup
            lastgroup = self.lastgroup
            groupdqflags = dqflags.group

            if firstgroup is not None or lastgroup is not None:
                set_groupdq(firstgroup, lastgroup, ngroups, result.groupdq, groupdqflags)

            # Before the ramp_fit() call, copy the input model ("_W" for weighting)
            # for later reconstruction of the fitting array tuples.
            input_model_w = result.copy()

            # Run ramp_fit(), ignoring all DO_NOT_USE groups, and return the
            # ramp fitting arrays for the ImageModel, the CubeModel, and the
            # RampFitOutputModel.
            image_info, integ_info, opt_info, gls_opt_model = ramp_fit.ramp_fit(
                result,
                buffsize,
                self.save_opt,
                readnoise_2d,
                gain_2d,
                self.algorithm,
                self.weighting,
                max_cores,
                dqflags.pixel,
                suppress_one_group=self.suppress_one_group,
            )

            # Create a gdq to modify if there are charge_migrated groups
            if self.algorithm == "OLS":
                gdq = input_model_w.groupdq.copy()

                dnu = dqflags.group["DO_NOT_USE"]
                chg = dqflags.group["CHARGELOSS"]
                sat = dqflags.group["SATURATED"]
                # Locate groups where that are flagged with CHARGELOSS
                wh_chargeloss = np.where(np.bitwise_and(gdq.astype(np.uint32), chg))

                if len(wh_chargeloss[0]) > 0:
                    # Unflag groups flagged as both CHARGELOSS and DO_NOT_USE
                    gdq[wh_chargeloss] -= dnu + chg

                    # Flag SATURATED groups as DO_NOT_USE for later segment determination
                    where_sat = np.where(np.bitwise_and(gdq, sat))
                    gdq[where_sat] = np.bitwise_or(gdq[where_sat], dnu)

                    # Get group_time for readnoise variance calculation
                    group_time = result.meta.exposure.group_time

                    # Using the modified GROUPDQ array, create new readnoise variance arrays
                    image_var_rn, integ_var_rn, opt_var_rn = compute_rn_variances(
                        gdq, readnoise_2d, gain_2d, group_time
                    )

                    # Create new ramp fitting array tuples, by inserting the new
                    # readnoise variances into copies of the original ramp fitting
                    # tuples.
                    image_info_new, integ_info_new = None, None
                    ch_int, ch_grp, ch_row, ch_col = wh_chargeloss
                    if image_info is not None and image_var_rn is not None:
                        rnoise = image_info[3]
                        rnoise[ch_row, ch_col] = image_var_rn[ch_row, ch_col]
                        image_info_new = (
                            image_info[0],
                            image_info[1],
                            image_info[2],
                            rnoise,
                            image_info[4],
                        )

                    if integ_info is not None and integ_var_rn is not None:
                        rnoise = integ_info[3]
                        rnoise[ch_int, ch_row, ch_col] = integ_var_rn[ch_int, ch_row, ch_col]
                        integ_info_new = (
                            integ_info[0],
                            integ_info[1],
                            integ_info[2],
                            rnoise,
                            integ_info[4],
                        )

                    image_info = image_info_new
                    integ_info = integ_info_new

                    opt_info_new = None
                    if opt_info is not None and opt_var_rn is not None:
                        opt_info_new = (
                            opt_info[0],
                            opt_info[1],
                            opt_info[2],
                            opt_var_rn,
                            opt_info[4],
                            opt_info[5],
                            opt_info[6],
                            opt_info[7],
                            opt_info[8],
                        )

                    opt_info = opt_info_new

        # Save the OLS optional fit product, if it exists.
        if opt_info is not None:
            opt_model = create_optional_results_model(result, opt_info)
            self.save_model(opt_model, "fitopt", output_file=self.opt_name)
        """
        # GLS removed from code, since it's not implemented right now.
        # Save the GLS optional fit product, if it exists
        if gls_opt_model is not None:
            self.save_model(
                gls_opt_model, 'fitoptgls', output_file=self.opt_name
            )
        """

        out_model, int_model = None, None
        # Create models from possibly updated info
        if image_info is not None and integ_info is not None:
            out_model = create_image_model(result, image_info)
            out_model.meta.bunit_data = "DN/s"
            out_model.meta.bunit_err = "DN/s"
            out_model.meta.cal_step.ramp_fit = "COMPLETE"
            if (result.meta.exposure.type in ["NRS_IFU", "MIR_MRS"]) or (
                result.meta.exposure.type in ["NRS_AUTOWAVE", "NRS_LAMP"]
                and result.meta.instrument.lamp_mode == "IFU"
            ):
                out_model = datamodels.IFUImageModel(out_model)

            int_model = create_integration_model(result, integ_info, int_times)
            int_model.meta.bunit_data = "DN/s"
            int_model.meta.bunit_err = "DN/s"
            int_model.meta.cal_step.ramp_fit = "COMPLETE"

        # Cleanup
        del result

        return out_model, int_model


def set_groupdq(firstgroup, lastgroup, ngroups, groupdq, groupdqflags):
    """
    Set the groupdq flags based on the values of firstgroup, lastgroup.

    Parameters
    ----------
    firstgroup : int or None
        The first group to be used in the ramp fitting
    lastgroup : int or None
        The last group to be used in the ramp fitting
    ngroups : int
        The number of groups in the ramp
    groupdq : ndarray
        The groupdq array to be modified in place
    groupdqflags : dict
        The groupdq flags dict
    """
    dnu = groupdqflags["DO_NOT_USE"]
    if firstgroup is None:
        firstgroup = 0

    if lastgroup is None:
        lastgroup = ngroups - 1

    if firstgroup < 0:
        log.warning("first group < 0, reset to 0")
        firstgroup = 0

    if lastgroup >= ngroups:
        log.warning(f"Last group number >= #groups ({ngroups}), reset to {ngroups - 1}")

    if firstgroup >= lastgroup:
        log.warning(f"firstgroup ({firstgroup}) cannot be >= lastgroup ({lastgroup})")
        log.warning("Group selectors ignored")
        firstgroup = 0
        lastgroup = ngroups - 1

    if firstgroup > 0:
        groupdq[:, :firstgroup] = np.bitwise_or(groupdq[:, :firstgroup], dnu)

    if lastgroup < (ngroups - 1):
        groupdq[:, (lastgroup + 1) :] = np.bitwise_or(groupdq[:, (lastgroup + 1) :], dnu)
