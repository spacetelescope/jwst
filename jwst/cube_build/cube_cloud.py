""" Map the detector pixels to the cube coordinate system.
This is where the weight functions are used.
"""
#from numba import jit
import numpy as np
import logging
from . import cube_overlap

#log = logging.getLogger('numba')
log.setLevel(logging.WARNING)


def find_closest_wave(iw,w,
                      wavelength_table,
                      rois_table,
                      roiw_table,
                      softrad_table,
                      weight_power_table,
                      scalerad_table,
                      rois_det,
                      roiw_det,
                      softrad_det,
                      weight_det,
                      scalerad_det):

    """ Given a specific wavelength, find the closest value in the wavelength_table

    """
    ifound = (np.abs(wavelength_table - w)).argmin()
    rois_det[iw] = rois_table[ifound]
    roiw_det[iw] = roiw_table[ifound]
    softrad_det[iw] = softrad_table[ifound]
    weight_det[iw] = weight_power_table[ifound]
    scalerad_det[iw] = scalerad_table[ifound]



