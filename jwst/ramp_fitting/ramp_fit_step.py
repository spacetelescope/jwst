#! /usr/bin/env python
import numpy as np

from stcal.ramp_fitting import ramp_fit

from stcal.ramp_fitting.likely_fit import LIKELY_MIN_NGROUPS


from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags


from jwst.stpipe import Step

from jwst.lib import reffile_utils

import logging

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


class RampFitStep(Step):
    """Fit line to determine the value of mean rate counts vs. time."""

    class_alias = "ramp_fit"

    spec = """
        algorithm = option('OLS_C', 'LIKELY', default='OLS_C')
        int_name = string(default='')
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        suppress_one_group = boolean(default=True)  # Suppress saturated ramps with good 0th group
        firstgroup = integer(default=None)   # Ignore groups before this one (zero indexed)
        lastgroup = integer(default=None)   # Ignore groups after this one (zero indexed)
        maximum_cores = string(default='1') # cores for multiprocessing. Can be an integer, 'half', 'quarter', or 'all'
    """  # noqa: E501

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

            int_times = result.int_times

            # Set the DO_NOT_USE bit in the groupdq values for groups before firstgroup
            # and groups after lastgroup
            firstgroup = self.firstgroup
            lastgroup = self.lastgroup
            groupdqflags = dqflags.group

            if firstgroup is not None or lastgroup is not None:
                set_groupdq(firstgroup, lastgroup, ngroups, result.groupdq, groupdqflags)

            # Run ramp_fit(), ignoring all DO_NOT_USE groups, and return the
            # ramp fitting arrays for the ImageModel, the CubeModel, and the
            # RampFitOutputModel.
            image_info, integ_info, opt_info = ramp_fit.ramp_fit(
                result,
                self.save_opt,
                readnoise_2d,
                gain_2d,
                self.algorithm,
                self.weighting,
                max_cores,
                dqflags.pixel,
                suppress_one_group=self.suppress_one_group,
            )

        # Save the OLS_C optional fit product, if it exists.
        if opt_info is not None:
            opt_model = create_optional_results_model(result, opt_info)
            self.save_model(opt_model, "fitopt", output_file=self.opt_name)

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
