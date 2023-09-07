import logging
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from gwcs.wcstools import grid_from_bounding_box

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Set of routines to find a  correction for a time-wavelength
# photometric response that is particularly significant at long wavelengths.


def expfunc_bounded(x, a, b, c, x0):
    """
    Short Summary
    --------------
    Time dependent flux loss function.

    The model parameters are a,b,c, x0. x0 is the reference day from which the
    time dependent parameters were derived. Each parameter is an array for
    different wavelength regions. This function will return an array of
    corrections at the given x day (MJD).

    Parameters
    ----------
    x : int
        Modified Julian Day (EXPMID) of the observation
    a : numpy float array
        Model parameter for time dependence
    b : numpy float array
        Model parameter for time dependence
    c : numpy float array
        Model parameter for time dependence
    x0 : int
        Reference day from which the time dependence
        coefficients were derived.
        (x-x0)/100 puts time in relevant units of 100 days

    Returns
    -------
    The time-dependent correction flux loss array.
    """

    temp = a * np.exp(-b * (x - x0) / 100) + c
    temp[temp > 1] = 1
    return temp


def get_correction_function(side, timecoeff, mid_time):
    """
    Short Summary
    -------------
    Constructing the time and wavelength dependent function

    The time dependence of the flux loss is described by the
    exponential function contained in expfunc_bounded.
    The MIRI MRS flux loss is time and wavelength dependent. The flux
    loss is much larger at the long wavelengths of the MRS.

    Parameters
    ----------
    side : str
        Either 'left' or 'right'. The MRS contains 2 channels in every
        exposure. The time-wavelength dependent correction is different
        for every MRS band.
    timecoeff : dictionary
        A dictionay holding the correction factors for the time/wavelength
        correction.  binwave, acoeff, bcoeff, ccoeff, x0 are the parameters
        for an MRS time-wavelength photom correction.
    mid-time : int
        Modified Julian day

    Returns
    -------
    Time-wavelength dependent photom loss correction

    """

    binwave = timecoeff[side]['binwave']
    a = timecoeff[side]['acoeff']
    b = timecoeff[side]['bcoeff']
    c = timecoeff[side]['ccoeff']
    x0 = timecoeff[side]['x0']

    # time dependent function. Time is set by mid_time (modified julian day)
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
    input : JWST IFUImageModel
        Input science data model to be corrected.

    detector : str
        MRS detector working on

    ftab : stdatamodels.jwst.datamodels.MirMrsPhotomModel
        MRS Photom reference file

    mid_time : float
        Exposure mid time in MJD

    Returns
    -------
    result : numpy array
        An array of corrections to apply to data.
    """

    # Read in the time-wavelength dependent coefficients
    # for each channel from the MRS photom reference file
    table_ch = {}
    table_ch['ch1'] = ftab.timecoeff_ch1
    table_ch['ch2'] = ftab.timecoeff_ch2
    table_ch['ch3'] = ftab.timecoeff_ch3
    table_ch['ch4'] = ftab.timecoeff_ch4

    # Each MIRI MRS exposure has 2 channels
    # We want to correct the channels seperately
    # If we have MIRIFUSHORT data then channel 1 is on the left
    # and channel 2 is on the right. If we have MIRIFULONG  data
    # then channel 4 is on the left and channel 3 is on right.

    timecoeff = {}
    timecoeff['left'] = {}
    timecoeff['right'] = {}
    left = 'ch1'
    right = 'ch2'
    timecoeff['left']['xstart'] = 0
    timecoeff['left']['xend'] = 512
    timecoeff['right']['xstart'] = 512
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

    ysize, xsize = input.data.shape[-2], input.data.shape[-1]
    # the correction is time and wavelength dependent. Pull out the
    # wavelength of the data
    x, y = grid_from_bounding_box(input.meta.wcs.bounding_box)
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
    result1[waveim == 0] = 0.0

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
    result2[waveim == 0] = 0.0

    # Combine the bands to have 1 correction for the entire detector
    result = np.zeros((ysize, xsize))
    result[:, l_xstart: l_xend] = result1
    result[:, r_xstart: r_xend] = result2
    indx = np.where(result == 0)
    result[indx] = 1.0
    result[np.isnan(result)] = 1.0

    return result
