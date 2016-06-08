#
#  Module for  subtracting lastframe correction from MIRI science data sets
#

import numpy as np
import logging
from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, lastframe_model):
    """
    Short Summary
    -------------
    Subtracts lastframe correction from science arrays, combines
    error arrays in quadrature, and updates data quality array based on 
    DQ flags in the lastframe arrays.

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    lastframe_model: lastframe model object
        lastframe data

    sci_ngroups:int
        number of groups in input data

    sci_nints: int
        number of integrations in input data

    Returns
    -------
    output_model: data model object
        lastframe-subtracted science data

    """

    # Save some data params for easy use later
    sci_nints   = input_model.data.shape[0]
    sci_ngroups = input_model.data.shape[1]

    #lastframe_nints = lastframe_model.data.shape[0]

    log.debug("LastFrame Sub using: nints=%d " %
          (sci_nints))



    # Create output as a copy of the input science data model
    output = input_model.copy()

    # combine the science and lastframe DQ arrays
    output.pixeldq = np.bitwise_or (input_model.pixeldq, lastframe_model.dq) 


    output.data[:,sci_ngroups-1] -= lastframe_model.data
           
       # combine the ERR arrays in quadrature
       # NOTE: currently stubbed out until ERR handling is decided
       # output.err[i,j] = np.sqrt( 
       # output.err[i,j]**2 + lastframe.err[j]**2 )
            
    return output


