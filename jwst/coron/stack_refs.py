"""
    Stack individual coronagraphic PSF reference images into a cube model.

:Authors: Howard Bushouse

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_
"""

import numpy as np
from jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def make_cube(input_table):
    """
    make_cube: Stack individual coronagraphic PSF reference images into a
    jwst_lib.CubeModel, for use in LOCI/KLIP processing.
    """

    # Get the number of input images
    num_refs = len(input_table.input_filenames)

    # Get the size of the first input image
    img = models.ImageModel(input_table.input_filenames[0])
    nrows, ncols = img.data.shape
    img.close()

    # Create an empty cube array of the appropriate dimensions
    cube = np.zeros((num_refs, nrows, ncols), dtype=np.float32)

    # Loop over the input images, copying only the science data array
    # into the cube
    for i, name in enumerate(input_table.input_filenames):
        log.info('Adding member %s', name)
        img = models.ImageModel(name)
        cube[i] = img.data
        img.close()

    # Create the ouput Cube model
    output_model = models.CubeModel(data=cube)

    return output_model

