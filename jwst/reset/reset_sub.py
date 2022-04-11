#
#  Module for  subtracting reset correction from  science data sets
#
import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, reset_model):
    """
    Short Summary
    -------------
    Subtracts reset correction from science arrays, combines
    error arrays in quadrature, and updates data quality array based on
    DQ flags in the reset arrays. When applying correction, if integration
    # of input data > reset_nints, use  correction of reset_nints to
    correct data. When apply correction, if group # of input data >
    reset_ngroups, use correction of reset_ngroups (which = 0).

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    reset_model: reset model object
        reset data


    Returns
    -------
    output_model: data model object
        reset-subtracted science data

    """

    # Save some data params for easy use later
    sci_nints = input_model.meta.exposure.nints      # num ints in input data
    sci_ngroups = input_model.meta.exposure.ngroups  # num groups in input data
    sci_integration_start = input_model.meta.exposure.integration_start
    sci_integration_end = input_model.meta.exposure.integration_end
    istart = 0
    iend = sci_nints

    if sci_integration_start is not None:
        istart = sci_integration_start - 1
    if sci_integration_end is not None:
        iend = sci_integration_end

    reset_nints = reset_model.meta.exposure.nints      # num of int in ref file
    reset_ngroups = reset_model.meta.exposure.ngroups  # num of groups

    # Replace NaN's in the reset with zeros(should not happen but just in case)
    reset_model.data[np.isnan(reset_model.data)] = 0.0
    log.debug("Reset Sub using: nints = {}, ngroups = {}".format(sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()

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
        output.pixeldq = np.bitwise_or(input_model.pixeldq, reset_model.dq)

        # we are only correcting the first reset_ngroups
        for j in range(igroup):
            output.data[i - istart, j] -= reset_model.data[ir, j]

            # combine the ERR arrays in quadrature
            # NOTE: currently stubbed out until ERR handling is decided
            # output.err[i,j] = np.sqrt(
            #           output.err[i,j]**2 + reset.err[j]**2)

    return output
