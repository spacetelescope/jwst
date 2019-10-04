from ..datamodels import dqflags
from ..lib import reffile_utils
from .linearity_func import apply_linearity_func

import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, lin_model):
    """
    Short Summary
    -------------
    Execute all tasks for linearity correction

    Parameters
    ----------
    input_model: data model object
        input science data model to be corrected

    lin_model: linearity model object
        linearity reference file data model

    Returns
    -------
    output_model: data model object
        linearity corrected data

    """
    # Create the output model as a copy of the input
    output_model = input_model.copy()

    # Propagate the DQ flags from the linearity ref data into the 2D science DQ
    propagate_dq_info(output_model, lin_model)

    # Apply the linearity correction coeffs to the science data
    apply_linearity(output_model, lin_model)

    return output_model


def propagate_dq_info(input, linearity_ref_model):
    """
    Short Summary
    -------------
    Propagate the data quality information from the linearity reference
    file DQ array to the 2D DQ array of the science data.

    Parameters
    ----------
    input: data model object
        science data model to be corrected in place

    linearity_ref_model: linearity model object
        linearity reference data

    Returns
    -------

    """

    # Retrieve 2D DQ subarray from linearity reference file, if necessary
    if reffile_utils.ref_matches_sci(input, linearity_ref_model):
        lindq = linearity_ref_model.dq
    else:
        log.info('Extracting linearity subarray to match science data')
        sub_lin_model = reffile_utils.get_subarray_model(input, linearity_ref_model)
        lindq = sub_lin_model.dq.copy()
        sub_lin_model.close()

    # Combine the DQ arrays using bitwise_or
    input.pixeldq = np.bitwise_or(input.pixeldq, lindq)


def apply_linearity(input, linearity_ref_model):
    """
    Short Summary
    -------------

    Apply the linearity correction to the pixels in the science ramp data
    that have not been flagged as saturated by the saturation step.

    Parameters
    ----------
    input: data model object
        The input science data to be corrected

    linearity_ref_model: linearity reference file model object

    Returns
    -------
    """

    ramp = input.data
    dq = input.groupdq

    # If the input data does not have an expanded DQ array, create one
    if len(dq) == 0:
        dq = (ramp * 0).astype(np.uint32)

    # Check for subarray mode
    if reffile_utils.ref_matches_sci(input, linearity_ref_model):
        lin_coeffs = linearity_ref_model.coeffs
        lin_dq = linearity_ref_model.dq
    else:
        sub_lin_model = reffile_utils.get_subarray_model(input, linearity_ref_model)
        lin_coeffs = sub_lin_model.coeffs.copy()
        lin_dq = sub_lin_model.dq.copy()

    # Check for NO_LIN_CORR flags in the DQ extension of the ref file
    lin_coeffs = correct_for_flag(lin_coeffs, lin_dq)

    # Check for NaNs in the COEFFS extension of the ref file
    lin_coeffs = correct_for_NaN(lin_coeffs, input)

    # Get the DQ bit value that represents saturation
    sat_val = dqflags.group['SATURATED']

    # Apply the correction function
    input.data = apply_linearity_func(ramp, dq, lin_coeffs, sat_val)


def correct_for_NaN(lin_coeffs, input):
    """
    Short Summary
    -------------
    Check for NaNs in the COEFFS extension of the ref file in case there are
    pixels that should have been (but were not) flagged there as NO_LIN_CORR
    (linearity correction not determined for pixel). For such pixels, update the
    coefficients so that there is effectively no correction, and flag their
    pixeldq values in place as NO_LIN_CORR in the step output.

    Parameters
    ----------
    lin_coeffs: 3D array
        array of correction coefficients in reference file

    input: data model object
        science data model to be corrected in place

    Returns
    -------
    lin_coeffs: 3D array
        updated array of correction coefficients in reference file
    """

    wh_nan = np.where(np.isnan(lin_coeffs))
    znan, ynan, xnan = wh_nan[0], wh_nan[1], wh_nan[2]
    num_nan = 0

    nan_array = np.zeros((lin_coeffs.shape[1], lin_coeffs.shape[2]),
                         dtype=np.uint32)

    # If there are NaNs as the correction coefficients, update those
    # coefficients so that those SCI values will be unchanged.
    if len(znan) > 0:
        ben_cor = ben_coeffs(lin_coeffs) # get benign coefficients
        num_nan = len(znan)

        for ii in range(num_nan):
            lin_coeffs[:, ynan[ii], xnan[ii]] = ben_cor
            nan_array[ynan[ii], xnan[ii]] = dqflags.pixel['NO_LIN_CORR']

        # Include these pixels in the output pixeldq
        input.pixeldq = np.bitwise_or(input.pixeldq, nan_array)

        log.debug("Unflagged pixels having coefficients set to NaN were"
                  " detected in the ref file; for those affected pixels"
                  " no linearity correction will be applied.")

    return lin_coeffs


def correct_for_flag(lin_coeffs, lin_dq):
    """
    Short Summary
    -------------
    Check for pixels that are flagged as NO_LIN_CORR
    ('No linearity correction available') in the DQ extension of the ref data.
    For such pixels, update the coefficients so that there is effectively no
    correction. Because these are already flagged in the ref file, they will
    also be flagged in the output dq.

    Parameters
    ----------
    lin_coeffs: 3D array
        array of correction coefficients in reference file

    lin_dq: 2D array
        array of data quality flags in reference file

    Returns
    -------
    lin_coeffs: 3D array
        updated array of correction coefficients in reference file
    """

    wh_flag = np.bitwise_and(lin_dq, dqflags.pixel['NO_LIN_CORR'])
    num_flag = len(np.where(wh_flag > 0)[0])

    wh_lin = np.where(wh_flag == dqflags.pixel['NO_LIN_CORR'])
    yf, xf = wh_lin[0], wh_lin[1]

    # If there are pixels flagged as 'NO_LIN_CORR', update the corresponding
    #     coefficients so that those SCI values will be unchanged.
    if (num_flag > 0):
        ben_cor = ben_coeffs(lin_coeffs) # get benign coefficients

        for ii in range(num_flag):
            lin_coeffs[:, yf[ii], xf[ii]] = ben_cor

        log.debug("Pixels were flagged in the DQ of the reference file as"
                  " NO_LIN_CORR ('Linearity correction not available'); for"
                  " those affected pixels no linearity correction will be"
                  " applied.")

    return lin_coeffs


def ben_coeffs(lin_coeffs):
    """
    Short Summary
    -------------
    For pixels having at least 1 NaN coefficient, reset the coefficients to be
    benign, which will effectively leave the SCI values unaffected.

    Parameters
    ----------
    lin_coeffs: 3D array
        array of correction coefficients in reference file

    Returns
    -------
    ben_cor: 1D array
        benign coefficients - all ben_cor[:] = 0.0 except ben_cor[1] = 1.0
    """
    ben_cor = np.zeros(lin_coeffs.shape[0])
    ben_cor[1] = 1.0

    return ben_cor
