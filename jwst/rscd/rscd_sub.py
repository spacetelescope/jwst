#
#  Module for applying the RSCD correction to science data
#

import numpy as np
import logging
from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, rscd_model):
    """
    Short Summary
    -------------
    Applies rscd correction to science arrays, xxx

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    rscd_model: rscd model object
        rscd data

    Returns
    -------
    output_model: data model object
        RSCD-corrected science data

    """

    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]       # number of integrations

    sci_ngroups = input_model.data.shape[1]     # number of groups
    # could also grab this information from input_model.meta.exposure.nints (ngroups)


    log.debug("RSCD correction using: nints=%d, ngroups=%d"  %
          (sci_nints, sci_ngroups))

    # Create output as a copy of the input science data model
    output = input_model.copy()

    # xxx Read the table in the reference file to get the coefficients.

    # loop over all integrations except the first

    for i in range(1, sci_nints):

       # loop over groups in input science data:
       for j in range(sci_ngroups):
           # Apply the correction
           output.data[i, j] += xxx

    return output
