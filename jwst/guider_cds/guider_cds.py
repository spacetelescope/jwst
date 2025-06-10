import warnings

import numpy as np
import logging

from stdatamodels.jwst import datamodels

from jwst.lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def guider_cds(model, gain_model, readnoise_model):
    """
    Calculate the count rate for each pixel in all integrations.

    For each integration in a given FGS guider dataset whose mode is ACQ1,
    ACQ2, or TRACK, the count rate is the last group minus the first group,
    divided by the effective integration time.  If the mode is ID, the last
    group minus the first group is calculated for both integrations; the count
    rate is then given by the minimum of these two values for each pixel,
    divided by the group time.  For the FINEGUIDE mode, the count rate is the
    average of the last 4 groups minus the average of the first 4 groups,
    divided by the group time.  The variances are calculated using gain and
    read noise values obtained from reference files if possible, otherwise,
    default representative constant values are used. From the variances the ERR
    array is calculated.

    Parameters
    ----------
    model : `datamodels.GuiderRawModel`
        Input data model
    gain_model : `datamodels.GainModel`
        Gain for all pixels
    readnoise_model : `datamodels.ReadnoiseModel`
        Readnoise for all pixels

    Returns
    -------
    new_model : 'datamodels.GuiderCalModel'
        Output data model
    """
    # get needed sizes and shapes
    grp_time = model.meta.exposure.group_time
    exp_type = model.meta.exposure.type
    n_int = model.data.shape[0]
    imshape = (model.data.shape[2], model.data.shape[3])

    # If exp_type is 'FGS_ID', force output to have single slice, and use
    # default GAIN and READNOISE values by setting their ref models to None
    if exp_type[:6] == "FGS_ID":
        new_model = datamodels.GuiderCalModel((1,) + imshape)
    else:
        new_model = datamodels.GuiderCalModel()

    # get gain and readnoise arrays to calculate ERR array
    gain_arr = get_ref_arr(model, gain_model)
    readnoise_arr = get_ref_arr(model, readnoise_model)

    # set up output data arrays
    slope_int_cube = np.zeros((n_int,) + imshape, dtype=np.float32)
    if exp_type[:6] == "FGS_ID":
        var_rn = np.zeros((1,) + imshape, dtype=np.float32)
        var_pn = np.zeros((1,) + imshape, dtype=np.float32)
    else:
        var_rn = np.zeros((n_int,) + imshape, dtype=np.float32)
        var_pn = np.zeros((n_int,) + imshape, dtype=np.float32)
    new_model.dq = model.dq

    # loop over data integrations
    for num_int in range(0, n_int):
        data_sect = model.data[num_int, :, :, :]
        if exp_type == "FGS_FINEGUIDE":
            first_4 = data_sect[:4, :, :].mean(axis=0)
            last_4 = data_sect[-4:, :, :].mean(axis=0)
            slope_int_cube[num_int, :, :] = last_4 - first_4

            var_rn[num_int, :, :] = 2 * (readnoise_arr / grp_time) ** 2
            var_pn[num_int, :, :] = slope_int_cube[num_int, :, :] / (gain_arr * grp_time)

        elif exp_type[:6] == "FGS_ID":
            grp_last = data_sect[1, :, :]
            grp_first = data_sect[0, :, :]

            if num_int == 0:
                diff_int0 = grp_last - grp_first
            if num_int == 1:
                diff_int1 = grp_last - grp_first

        else:  # ACQ1, ACQ2, or TRACK
            grp_last = data_sect[1, :, :]
            grp_first = data_sect[0, :, :]
            slope_int_cube[num_int, :, :] = grp_last - grp_first
            var_rn[num_int, :, :] = 2 * (readnoise_arr / grp_time) ** 2
            var_pn[num_int, :, :] = slope_int_cube[num_int, :, :] / (gain_arr * grp_time)

    if exp_type[:6] == "FGS_ID":
        new_model.data[0, :, :] = np.minimum(diff_int1, diff_int0) / grp_time
        var_rn[0, :, :] = 2 * (readnoise_arr / grp_time) ** 2

        # May be zeros in gain array - var will be NaN in these pixels.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "divide by zero", RuntimeWarning)
            warnings.filterwarnings("ignore", "invalid value", RuntimeWarning)
            var_pn[0, :, :] = np.minimum(diff_int1, diff_int0) / (gain_arr * grp_time)

    else:  # FINEGUIDE, ACQ1, ACQ2, or TRACK
        new_model.data = slope_int_cube / grp_time

    # set err to sqrt of sum of variances
    var_pn[var_pn < 0] = 0.0  # ensure variance is non-negative
    new_model.err = (var_rn + var_pn) ** 0.5

    # Add all table extensions to be carried over to output
    if len(model.planned_star_table):
        new_model.planned_star_table = model.planned_star_table
    if len(model.flight_star_table):
        new_model.flight_star_table = model.flight_star_table
    if len(model.pointing_table):
        new_model.pointing_table = model.pointing_table
    if len(model.centroid_table):
        new_model.centroid_table = model.centroid_table
    if len(model.track_sub_table):
        new_model.track_sub_table = model.track_sub_table

    # copy all meta data from input to output model
    new_model.update(model)

    # Update BUNIT to reflect count rate
    new_model.meta.bunit_data = "DN/s"

    return new_model


def get_ref_arr(model, reference_model):
    """
    Extract arrays in appropriate shape from the reference files.

    Parameters
    ----------
    model : `datamodels.GuiderRawModel`
        Input data model
    reference_model : `datamodels.ReferenceFileModel`
        Reference file containing the relevant array;
        typically either a GainModel or a ReadnoiseModel

    Returns
    -------
    ref_arr : ndarray, 2-D, float
        Values from the reference file
    """
    # extract subarray from reference file, if necessary
    if reffile_utils.ref_matches_sci(model, reference_model):
        ref_arr = reference_model.data
    else:
        try:
            ref_sub_model = reffile_utils.get_subarray_model(model, reference_model)
            log.info(
                f"Extracting subarray from reference {reference_model.meta.model_type} "
                f"to match science data."
            )
            ref_arr = ref_sub_model.data
            ref_sub_model.close()
        except ValueError as err:
            if "STACK" in model.meta.exposure.type:
                if model.data.shape[-2:] != reference_model.data.shape:
                    log.debug(
                        f"The {reference_model.meta.model_type} reference file array does not "
                        "match the shape of stacked FGS data;"
                        "applying the mean value to the data."
                    )
                    ref_arr = np.zeros(model.data.shape[-2:], dtype=np.float32) + np.mean(
                        reference_model.data
                    )
            else:
                raise ValueError from err

    return ref_arr
