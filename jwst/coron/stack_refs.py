"""
    Stack individual coronagraphic PSF reference images into a cube model.

:Authors: Howard Bushouse

"""

import numpy as np
from stdatamodels.jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_cube(input_models):
    """
    make_cube: Stack all of the integrations from multiple PSF
    reference exposures into a single CubeModel, for use in the
    coronagraphic alignment and PSF-subtraction steps.
    """

    # Get the number of input images
    num_refs = len(input_models)

    # Loop over all the inputs to find the total number of integrations
    nrows_ref, ncols_ref = input_models[0].shape[-2:]
    nints = 0
    for i in range(num_refs):
        nints += input_models[i].shape[0]
        nrows, ncols = input_models[i].shape[-2:]
        if nrows != nrows_ref or ncols != ncols_ref:
            raise ValueError('All PSF exposures must have the same x/y dimensions!')

    # Create empty output data arrays of the appropriate dimensions
    outdata = np.zeros((nints, nrows, ncols), dtype=np.float64)
    outerr = outdata.copy()
    outdq = np.zeros((nints, nrows, ncols), dtype=np.uint32)

    # Loop over the input images, copying the data arrays
    # into the output arrays
    nint = 0
    for i in range(num_refs):
        log.info(' Adding psf member %d to output stack', i + 1)
        for j in range(input_models[i].shape[0]):
            outdata[nint] = input_models[i].data[j]
            outerr[nint] = input_models[i].err[j]
            outdq[nint] = input_models[i].dq[j]
            nint += 1

    # Create the output Cube model
    output_model = datamodels.CubeModel(data=outdata, err=outerr, dq=outdq)
    output_model.update(input_models[0])  # copy input meta data to output

    return output_model
