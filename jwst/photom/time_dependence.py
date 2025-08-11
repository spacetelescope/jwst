import logging

import numpy as np
from gwcs.wcstools import grid_from_bounding_box
from scipy.interpolate import RegularGridInterpolator

log = logging.getLogger(__name__)

__all__ = [
    "linear_correction",
    "exponential_correction",
    "powerlaw_correction",
    "get_correction_table",
    "miri_mrs_time_correction",
]


def linear_correction(midtime, t0, lossperyear, bounded=False):
    """
    Linear correction to photometry value.

    Parameters
    ----------
    midtime : float
        Mid-point MJD of observation.
    t0 : ndarray of float
        Reference day in MJD.
    lossperyear : ndarray of float
        Fractional loss of throughput per year (linear slope).
    bounded : bool, optional
        If True, any correction value greater than 1 is set to 1.0.

    Returns
    -------
    correction : ndarray of float
        Multiplicative correction values corresponding to the input parameters.
    """
    correction = 1.0 - lossperyear * (midtime - t0) / 365.0
    if bounded:
        correction[correction > 1.0] = 1.0
    return correction


def exponential_correction(midtime, t0, amplitude, tau, const, bounded=False):
    """
    Exponential correction to photometry value.

    Parameters
    ----------
    midtime : float
        Mid-point MJD of observation.
    t0 : ndarray of float
        Reference day in MJD.
    amplitude : ndarray of float
        Exponential amplitude.
    tau : ndarray of float
        E-folding time constant.
    const : ndarray of float
        Long-term exponential asymptote.
    bounded : bool, optional
        If True, any correction value greater than 1 is set to 1.0.

    Returns
    -------
    correction : ndarray of float
        Multiplicative correction values corresponding to the input parameters.
    """
    correction = amplitude * np.exp(-(midtime - t0) / tau) + const
    if bounded:
        correction[correction > 1.0] = 1.0
    return correction


def powerlaw_correction(midtime, t0, year1value, tsoft, alpha, bounded=False):
    """
    Power law correction to photometry value.

    Parameters
    ----------
    midtime : float
        Mid-point MJD of observation.
    t0 : ndarray of float
        Reference day in MJD.
    year1value : ndarray of float
        Relative throughput 1 year after t0.
    tsoft : ndarray of float
        Softening parameter for power law decline.
    alpha : ndarray of float
        Power law loss coefficient.
    bounded : bool, optional
        If True, any correction value greater than 1 is set to 1.0.

    Returns
    -------
    correction : ndarray of float
        Multiplicative correction values corresponding to the input parameters.
    """
    norm = ((365.0 + tsoft) ** alpha) / year1value
    correction = ((midtime - t0 + tsoft) ** alpha) / norm
    if bounded:
        correction[correction > 1.0] = 1.0
    return correction


def get_correction_table(photom_model, midtime, bounded=False):
    """
    Get a time-dependence correction table from a PHOTOM reference file model.

    If a ``phot_table`` is present, the table returned matches the
    shape of the phot_table.  If it is not present, ``None`` is returned.

    Time correction parameters are expected to be present in
    ``timecoeff_linear``, ``timecoeff_exponential``, or ``timecoeff_powerlaw``
    attributes. If no time correction parameters are present
    in the model, all values in the output array are 1.0.  If correction
    parameters are present, the correction value is computed from
    the expected functional form and the input time.  If more than one
    correction is present, they are multiplied together.

    Parameters
    ----------
    photom_model : DataModel
        May be any photom datamodel. Expected attributes are phot_table,
        timecoeff_linear, timecoeff_exponential, and timecoeff_powerlaw.
    midtime : float
        Mid-point MJD of observation.
    bounded : bool, optional
        If True, any correction value greater than 1 is set to 1.0.

    Returns
    -------
    correction : ndarray of float or None
        None is returned if the input model had no phot_table. Otherwise,
        the correction array matches the length of the input phot_table.
    """
    if not hasattr(photom_model, "phot_table"):
        return None

    # Multiply together any existing corrections
    correction = np.full(photom_model.phot_table.shape, 1.0)
    if photom_model.hasattr("timecoeff_linear"):
        param = photom_model.timecoeff_linear
        correction *= linear_correction(midtime, param["t0"], param["lossperyear"], bounded=bounded)
    if photom_model.hasattr("timecoeff_exponential"):
        param = photom_model.timecoeff_exponential
        correction *= exponential_correction(
            midtime, param["t0"], param["amplitude"], param["tau"], param["const"], bounded=bounded
        )
    if photom_model.hasattr("timecoeff_powerlaw"):
        param = photom_model.timecoeff_powerlaw
        correction *= powerlaw_correction(
            midtime,
            param["t0"],
            param["year1value"],
            param["tsoft"],
            param["alpha"],
            bounded=bounded,
        )

    return correction


def _get_mrs_correction_function(side, timecoeff, mid_time):
    """
    Construct the time and wavelength dependent correction function.

    The MIRI MRS flux loss is time and wavelength dependent. The flux
    loss is much larger at the long wavelengths of the MRS.

    Parameters
    ----------
    side : str
        Either "left" or "right". The MRS contains 2 channels in every
        exposure. The time-wavelength dependent correction is different
        for every MRS band.
    timecoeff : dictionary
        A dictionary holding the correction factors for the time/wavelength
        correction.  binwave, acoeff, bcoeff, ccoeff, x0 are the parameters
        for an MRS time-wavelength photom correction.
    mid_time : float
        Modified Julian day

    Returns
    -------
    function
        Time-wavelength dependent photom loss correction
    """
    binwave = timecoeff[side]["binwave"]
    a = timecoeff[side]["acoeff"]
    b = timecoeff[side]["bcoeff"]
    c = timecoeff[side]["ccoeff"]
    x0 = timecoeff[side]["x0"]

    # Timescale is (100/b) days
    tau = 100 / b

    # Time dependent function. Time is set by mid_time (modified julian day)
    corgrid = exponential_correction(mid_time, x0, a, tau, c, bounded=True)

    # Now fold in the wavelength dependence
    func = RegularGridInterpolator([binwave], corgrid, bounds_error=False, fill_value=None)
    return func


def miri_mrs_time_correction(input_model, detector, ftab, mid_time):
    """
    Find the time and wavelength dependent photom correction for MRS data.

    Parameters
    ----------
    input_model : JWST IFUImageModel
        Input science data model to be corrected.
    detector : str
        MRS detector working on
    ftab : MirMrsPhotomModel
        MRS Photom reference file
    mid_time : float
        Exposure mid time in MJD

    Returns
    -------
    result : numpy.ndarray
        An array of corrections to apply to data.
    """
    # Read in the time-wavelength dependent coefficients
    # for each channel from the MRS photom reference file
    table_ch = {}
    table_ch["ch1"] = ftab.timecoeff_ch1
    table_ch["ch2"] = ftab.timecoeff_ch2
    table_ch["ch3"] = ftab.timecoeff_ch3
    table_ch["ch4"] = ftab.timecoeff_ch4
    # Each MIRI MRS exposure has 2 channels
    # We want to correct the channels separately
    # If we have MIRIFUSHORT data then channel 1 is on the left
    # and channel 2 is on the right. If we have MIRIFULONG  data
    # then channel 4 is on the left and channel 3 is on right.

    timecoeff = {}
    timecoeff["left"] = {}
    timecoeff["right"] = {}
    left = "ch1"
    right = "ch2"
    timecoeff["left"]["xstart"] = 0
    timecoeff["left"]["xend"] = 512
    timecoeff["right"]["xstart"] = 512
    timecoeff["right"]["xend"] = 1031
    if detector == "MIRIFULONG":
        left = "ch4"
        right = "ch3"
        # Check the reference file has the time dependent coefficients
        # check that table 1 wavelength bin is an array with values

    # Pull out the time coefficients for the detector we are working on
    timecoeff["left"]["binwave"] = table_ch[left]["binwave"]
    timecoeff["left"]["acoeff"] = table_ch[left]["acoeff"]
    timecoeff["left"]["bcoeff"] = table_ch[left]["bcoeff"]
    timecoeff["left"]["ccoeff"] = table_ch[left]["ccoeff"]
    timecoeff["left"]["x0"] = table_ch[left]["x0"]

    timecoeff["right"]["binwave"] = table_ch[right]["binwave"]
    timecoeff["right"]["acoeff"] = table_ch[right]["acoeff"]
    timecoeff["right"]["bcoeff"] = table_ch[right]["bcoeff"]
    timecoeff["right"]["ccoeff"] = table_ch[right]["ccoeff"]
    timecoeff["right"]["x0"] = table_ch[right]["x0"]

    ysize, xsize = input_model.data.shape[-2], input_model.data.shape[-1]
    # the correction is time and wavelength dependent. Pull out the
    # wavelength of the data
    x, y = grid_from_bounding_box(input_model.meta.wcs.bounding_box)
    _, _, wave = input_model.meta.wcs(x, y)

    # correction for left side: set up the pixels to extract this region
    side = "left"
    l_xstart = timecoeff[side]["xstart"]
    l_xend = timecoeff[side]["xend"]

    # Find the function based on time and wavelength
    func = _get_mrs_correction_function(side, timecoeff, mid_time)
    waveim = wave[:, l_xstart:l_xend]

    # Determine the correction based input wavelength
    result1_1d = func(waveim.ravel())
    result1 = np.reshape(result1_1d, waveim.shape)
    result1[waveim == 0] = 0.0

    # correction of the right side: set up the pixels to extract this region
    side = "right"
    r_xstart = timecoeff[side]["xstart"]
    r_xend = timecoeff[side]["xend"]

    # Find the function based on time and wavelength
    func = _get_mrs_correction_function(side, timecoeff, mid_time)
    waveim = wave[:, r_xstart:r_xend]

    # Determine the correction based on input wavelength
    result2_1d = func(waveim.ravel())
    result2 = np.reshape(result2_1d, waveim.shape)
    result2[waveim == 0] = 0.0

    # Combine the bands to have 1 correction for the entire detector
    result = np.zeros((ysize, xsize))
    result[:, l_xstart:l_xend] = result1
    result[:, r_xstart:r_xend] = result2
    indx = np.where(result == 0)
    result[indx] = 1.0
    result[np.isnan(result)] = 1.0

    return result
