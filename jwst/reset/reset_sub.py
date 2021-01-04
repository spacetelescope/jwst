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
    sci_nints = input_model.data.shape[0]           # number of integrations in input data

    sci_ngroups = input_model.data.shape[1]           # number of groups in input data
    # could also grab this information from input_model.meta.exposure.nints (ngroups)

    reset_nints = reset_model.data.shape[0]           # number of integrations in reset reference file
    reset_ngroups = reset_model.data.shape[1]         # number of groups

     # Replace NaN's in the reset with zeros (should not happen but just in case)
    reset_model.data[np.isnan(reset_model.data)] = 0.0
    log.debug("Reset Sub using: nints=%d, ngroups=%d" %
          (sci_nints, sci_ngroups))


    # Create output as a copy of the input science data model
    output = input_model.copy()

    # loop over all integrations
    for i in range(sci_nints):
        # check if integration # > reset_nints
        ir = i

        if i >= reset_nints:
            ir = reset_nints - 1

        # combine the science and reset DQ arrays
        output.pixeldq = np.bitwise_or(input_model.pixeldq, reset_model.dq)

        # loop over groups in input science data:
        for j in range(sci_ngroups):
            jr = j
            if j <= (reset_ngroups-1):
                # subtract the SCI arrays for the groups = < reset_ngroups
                output.data[i, j] -= reset_model.data[ir, jr]

            # combine the ERR arrays in quadrature
            # NOTE: currently stubbed out until ERR handling is decided
            #output.err[i,j] = np.sqrt(
            #           output.err[i,j]**2 + reset.err[j]**2)


    return output
