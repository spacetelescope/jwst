import numpy as np

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def subtract(model1, model2):
    """
    Short Summary
    -------------
    Subtract one data model from another, and include updated DQ
    in output.

    Parameters
    ----------
    model1: JWST data model
        input data model on which subtraction will be performed

    model2: JWST data model
        input data model that will be subtracted from the first model

    Returns
    -------
    output: JWST data model
        subtracted data model
    """

    # Create the output model as a copy of the first input
    output = model1.copy()

    # Subtract the SCI arrays
    output.data = model1.data - model2.data

    # Combine the ERR arrays in quadrature
    # NOTE: currently stubbed out until ERR handling is decided
    # output.err = np.sqrt(model1.err**2 + model2.err**2)

    # Combine the DQ flag arrays using bitwise OR
    output.dq = np.bitwise_or(model1.dq, model2.dq)

    # Return the subtracted model
    return output


def subtract_ints(sci_mod, bkg_mod):
    """
    Short Summary
    -------------
    Subtract averaged bkg model from sci model for rateints input, and
    include updated DQ in output.

    Parameters
    ----------
    sci_mod: JWST Cube data model for rateints (3D) SCI exposure
        Input data model on which subtraction will be performed. The 0th axis
        represents the integrations of the exposure.

    bkg_mod: JWST Image data model for the input BKG exposures.
        This BKG image has already been averaged over all rateint (3D)
        background exposures and all their integrations.

    Returns
    -------
    output: JWST Cube data model
        subtracted data model
    """
    # Create the output model as a copy of the SCI exposure
    output = sci_mod.copy()

    # Subtract the background average of data from each integration of the SCI arrays
    output.data = sci_mod.data - bkg_mod.data

    # bitwise_or combine DQ in BKG average with each SCI integration
    for i_nint in range(sci_mod.dq.shape[0]):
        output.dq[i_nint, :, :] = np.bitwise_or(output.dq[i_nint, :, :], bkg_mod.dq[:, :])

    # Combine the ERR arrays in quadrature
    # NOTE: currently stubbed out until ERR handling is decided
    # output.err = np.sqrt(sci_mod.err**2 + bkg_mod.err**2)
    # This is the same behavior as currently in code for the rate input

    # Return the subtracted model
    return output
