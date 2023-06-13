import logging
import numpy as np
from scipy.interpolate import RegularGridInterpolator

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def expfunc_bounded(x, a, b, c, x0):
    """ Function using time dependent coefficients """

    temp = a * np.exp(-b * (x - x0) / 100) + c
    temp[temp > 1] = 1
    return temp


def get_correction_function(side, timecoeff, mid_time):
    """ Find the time and wavelength dependent function """

    binwave = timecoeff[side]['binwave']
    a = timecoeff[side]['acoeff']
    b = timecoeff[side]['bcoeff']
    c = timecoeff[side]['ccoeff']
    x0 = timecoeff[side]['x0']

    # time depenpent function
    corgrid = expfunc_bounded(mid_time, a, b, c, x0)
    # now fold in the wavelength dependence
    func = RegularGridInterpolator([binwave], corgrid,
                                   bounds_error=False, fill_value=None)
    return func


def time_correction(input, detector, ftab, mid_time):
    """
    Find the time and wavelength dependent photom correction for MRS data.

    Parameters
    ----------
    input : JWST data model
        Input science data model to be corrected.

    detector : string
        MRS detector working on

    ftab : stdatamodels.jwst.datamodels.MirMrsPhotomModel
        MRS PHotom reference file

    mid_time: float
        Exposure mid time in MJD

    Returns
    -------
    result : numpy array
        An array of corrections to apply to data.
    """

    # Tables to be read in from MRS Photom reference file
    # Read in the time coefficients for each channel
    table_ch = {}
    table_ch['ch1'] = ftab.timecoeff_ch1
    table_ch['ch2'] = ftab.timecoeff_ch2
    table_ch['ch3'] = ftab.timecoeff_ch3
    table_ch['ch4'] = ftab.timecoeff_ch4

    # Each MIRI MRS exposure has 2 channels
    # We want to correct the channels seperately
    # If we have MIRIFUSHORT data then channel 1 is on the left and channel 2
    # is on the right. If we have MIRIFULONG  data then channel 4 is on the left
    # and channel 3 is on right.

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

    # Pull out the time coefficents for the detector we are working on
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
    # the correction is time and wavelength dependent. Pull out the
    # wavelength of the data
    _, _, wave = input.meta.wcs(x, y)

    # correction for left side: set up the pixels to extract this region
    side = 'left'
    l_xstart = timecoeff[side]['xstart']
    l_xend = timecoeff[side]['xend']

    # Find the function based on time and wavelength
    func = get_correction_function(side, timecoeff, mid_time)
    waveim = wave[:, l_xstart:l_xend]

    # Determine the correction based input wavelength
    result1_1d = func(waveim.ravel())
    result1 = np.reshape(result1_1d, waveim.shape)
    indx = np.where(waveim == 0)
    result1[indx] = 0.0

    # correction of the right side: set up the pixels to extract this region
    side = 'right'
    r_xstart = timecoeff[side]['xstart']
    r_xend = timecoeff[side]['xend']

    # Find the function based on time and wavelength
    func = get_correction_function(side, timecoeff, mid_time)
    waveim = wave[:, r_xstart: r_xend]

    # Determine the correction based on input wavelength
    result2_1d = func(waveim.ravel())
    result2 = np.reshape(result2_1d, waveim.shape)
    indx = np.where(waveim == 0)
    result2[indx] = 0.0

    # Combine the bands to have 1 correction for the entire detector
    result = np.zeros((ysize, xsize))
    result[:, l_xstart: l_xend] = result1
    result[:, r_xstart: r_xend] = result2
    indx = np.where(result == 0)
    result[indx] = 1.0
    indx = np.where(np.isnan(result))
    result[indx] = 1.0

    return result
