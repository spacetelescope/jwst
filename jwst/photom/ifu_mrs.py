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


def time_correction(input, detector, ftab, mid_time):

    table_ch = {}
    table_ch['ch1'] = ftab.timecoeff_ch1
    table_ch['ch2'] = ftab.timecoeff_ch2
    table_ch['ch3'] = ftab.timecoeff_ch3
    table_ch['ch4'] = ftab.timecoeff_ch4

    timecoeff = {}
    timecoeff['left'] = {}
    timecoeff['right'] = {}
    left = 'ch1'
    right = 'ch2'
    timecoeff['left']['xstart'] = 0
    timecoeff['left']['xend'] = 512
    timecoeff['right']['xstart'] = 513
    timecoeff['right']['xend'] = 1031
    if detector == 'MIRIFULONG':
        left = 'ch4'
        right = 'ch3'
        # Check the reference file has the time dependent coefficients
        # check that table 1 wavelength bin is an array with values

    timecoeff['left']['binwave'] = table_ch[left]['binwave']
    timecoeff['left']['acoeff'] = table_ch[left]['acoeff']
    timecoeff['left']['bcoeff'] = table_ch[left]['bcoeff']
    timecoeff['left']['ccoeff'] = table_ch[left]['ccoeff']
    timecoeff['left']['x0'] = table_ch[left]['x0']

    timecoeff['right']['binwave'] = table_ch[right]['binwave']
    timecoeff['right']['acoeff'] = table_ch[right]['acoeff']
    timecoeff['right']['bcoeff'] = table_ch[right]['bcoeff']
    timecoeff['right']['ccoeff'] = table_ch[right]['ccoeff']
    timecoeff['right']['x0'] = table_ch[right]['x0']

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
