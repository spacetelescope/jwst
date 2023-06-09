import logging
import functools
import warnings

import numpy as np
from astropy import units as u

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from scipy.interpolate import RegularGridInterpolator

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def expfunc_bounded(x, a, b, c, x0):
    temp = a * np.exp(-b * (x - x0) / 100) + c
    temp[temp > 1] = 1
    return temp


def get_correction_function(side, timecoeff, mid_time):
    binwave = timecoeff[side]['binwave']
    a = timecoeff[side]['acoeff']
    b = timecoeff[side]['bcoeff']
    c = timecoeff[side]['ccoeff']
    x0 = timecoeff[side]['x0']

    nbins = len(a)
    corgrid = expfunc_bounded(mid_time, a, b, c, x0)
    func = RegularGridInterpolator([binwave], corgrid,
                                   bounds_error=False, fill_value=None)
    return func


def time_correction(input, timecoeff, mid_time):

    ysize, xsize = input.data.shape
    y, x = np.mgrid[:ysize, :xsize]
    _, _, wave = input.meta.wcs(x, y)

    # correction for left side
    side = 'left'
    l_xstart = timecoeff[side]['xstart']
    l_xend = timecoeff[side]['xend']

    func = get_correction_function(side, timecoeff, mid_time)
    waveim = wave[:, l_xstart:l_xend]
    result1_1d = func(waveim.ravel())
    result1 = np.reshape(result1_1d, waveim.shape)
    indx = np.where(waveim == 0)
    result1[indx] = 0.0

    # Do the right side
    side = 'right'
    r_xstart = timecoeff[side]['xstart']
    r_xend = timecoeff[side]['xend']

    func = get_correction_function(side, timecoeff, mid_time)
    waveim = wave[:, r_xstart: r_xend]

    result2_1d = func(waveim.ravel())
    result2 = np.reshape(result2_1d, waveim.shape)
    indx = np.where(waveim == 0)
    result2[indx] = 0.0

    # Combine the bands
    result = np.zeros((ysize, xsize))
    result[:, l_xstart: l_xend] = result1
    result[:, r_xstart: r_xend] = result2
    indx = np.where(result == 0)
    result[indx] = 1.0
    indx = np.where(np.isnan(result))
    result[indx] = 1.0

    return result
