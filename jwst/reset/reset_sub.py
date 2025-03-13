#
#  Module for  subtracting reset correction from  science data sets
#
import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(output_model, reset_model):
    """
    Subtract the reset correction from science arrays.

    Subtracts the reset correction from science arrays and updates data quality array based
    on DQ flags in the reset arrays. The reset data model contains the number of initial
    frames in an integration to correct in the reset_ngroups parameter and the number of
    integrations for which a correction exists in the reset_nints parameter. Only the
    first reset_ngroups in each integration will have a reset correction applied. If the
    number of integrations in the science data is larger than reset_nints, then the reset correction
    applied is the correction contained  in the reset_nints integration.

    Parameters
    ----------
    output_model : `~jwst.datamodel.RampModel`
        Science data to be corrected

    reset_model : `~jwst.datamodel.ResetModel`
        Reset reference file model

    Returns
    -------
    output_model : `~jwst.datamodel.RampModel`
        Reset-subtracted science data
    """
    # Save some data params for easy use later
    sci_nints = output_model.meta.exposure.nints  # num ints in input data
    sci_ngroups = output_model.meta.exposure.ngroups  # num groups in input data
    sci_integration_start = output_model.meta.exposure.integration_start
    sci_integration_end = output_model.meta.exposure.integration_end
    istart = 0
    iend = sci_nints

    if sci_integration_start is not None:
        istart = sci_integration_start - 1
    if sci_integration_end is not None:
        iend = sci_integration_end

    reset_nints = reset_model.meta.exposure.nints  # num of int in ref file
    reset_ngroups = reset_model.meta.exposure.ngroups  # num of groups

    # Replace NaN's in the reset with zeros(should not happen but just in case)
    reset_model.data[np.isnan(reset_model.data)] = 0.0
    log.debug(f"Reset Sub using: nints = {sci_nints}, ngroups = {sci_ngroups}")

    # find out how many groups  we are correcting
    # the maximum number of groups to correct is reset_ngroups
    igroup = sci_ngroups
    if reset_ngroups < sci_ngroups:
        igroup = reset_ngroups

    # loop over all integrations in science image
    for i in range(istart, iend):
        # check if integration # > reset_nints
        ir = i
        if i >= reset_nints:
            ir = reset_nints - 1

        # combine the science and reset DQ arrays
        output_model.pixeldq = np.bitwise_or(output_model.pixeldq, reset_model.dq)

        # we are only correcting the first reset_ngroups
        for j in range(igroup):
            output_model.data[i - istart, j] -= reset_model.data[ir, j]

            # combine the ERR arrays in quadrature
            # NOTE: currently stubbed out until ERR handling is decided
            # output.err[i,j] = np.sqrt(
            #           output.err[i,j]**2 + reset.err[j]**2)

    return output_model
