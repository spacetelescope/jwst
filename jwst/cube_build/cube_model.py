# Routines used for building cubes
from __future__ import absolute_import, print_function

import sys
import time
import numpy as np
import math
import json

from astropy.io import fits

from .. import datamodels
from ..assign_wcs import pointing
from . import cube


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


#********************************************************************************

def SetUpIFUCube(Cube):

#********************************************************************************
    """
    Short Summary
    -------------
    Write the IFU cube to fits file

    Parameters
    ----------
    Cube: holds meta data of cube
    spaxel: list of spaxels in cube


    Returns
    -------
    no return = writes file

    """

    data = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))
    idata = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))

    dq_cube = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))
    err_cube = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))

    IFUCube = datamodels.IFUCubeModel(data=data, dq=dq_cube, err=err_cube, weightmap=idata)

    IFUCube.meta.filename = Cube.output_name
    IFUCube.meta.wcsinfo.crval1 = Cube.Crval1
    IFUCube.meta.wcsinfo.crval2 = Cube.Crval2
    IFUCube.meta.wcsinfo.crval3 = Cube.Crval3
    IFUCube.meta.wcsinfo.crpix1 = Cube.Crpix1
    IFUCube.meta.wcsinfo.crpix2 = Cube.Crpix2
    IFUCube.meta.wcsinfo.crpix3 = Cube.Crpix3
    IFUCube.meta.wcsinfo.cdelt1 = Cube.Cdelt1/3600.0
    IFUCube.meta.wcsinfo.cdelt2 = Cube.Cdelt2/3600.0
    IFUCube.meta.wcsinfo.cdelt3 = Cube.Cdelt3/3600.0

    IFUCube.meta.wcsinfo.ctype1 = 'RA---TAN'
    IFUCube.meta.wcsinfo.ctype2 = 'DEC--TAN'
    IFUCube.meta.wcsinfo.ctype3 = 'WAVE'

    IFUCube.meta.wcsinfo.cunit1 = 'deg'
    IFUCube.meta.wcsinfo.cunit2 = 'deg'
    IFUCube.meta.wcsinfo.cunit3 = 'um'

    IFUCube.meta.wcsinfo.wcsaxes = 3

    IFUCube.meta.flux_extension = 'SCI'
    IFUCube.meta.error_extension = 'ERR'
    IFUCube.meta.dq_extension = 'DQ'
    IFUCube.meta.weightmap = 'WMAP'
    IFUCube.meta.data_model_type = 'IFUCubeModel'
    IFUCube.error_type = 'ERR'


    wcsobj = pointing.create_fitswcs(IFUCube)
    IFUCube.meta.wcs = wcsobj

    return IFUCube


#********************************************************************************
def UpdateIFUCube(self, Cube,IFUCube, spaxel):

#********************************************************************************
    """
    Short Summary
    -------------
    Write the IFU cube to fits file

    Parameters
    ----------
    Cube: holds meta data of cube
    spaxel: list of spaxels in cube


    Returns
    -------
    no return = writes file

    """
    #pull out data into array


    icube = 0
    for z in range(Cube.naxis3):
        for y in range(Cube.naxis2):
            for x in range(Cube.naxis1):
                IFUCube.data[z, y, x] = spaxel[icube].flux
                IFUCube.weightmap[z, y, x] = len(spaxel[icube].ipointcloud)
                icube = icube + 1



    IFUCube.meta.filename = Cube.output_name   

    IFUCube.save(IFUCube.meta.filename)
    IFUCube.close()


    log.info('Wrote %s', IFUCube.meta.filename)
    return IFUCube

#********************************************************************************



