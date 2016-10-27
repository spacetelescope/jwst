from __future__ import division

#
#  Module for calculating pathloss correction for science data sets
#

import numpy as np
import logging
from .. import datamodels
from jwst.assign_wcs import nirspec

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def getCenter(exp_type, input):
    """
    Get the center of the target in the aperture.
    (0.0, 0.0) is the aperture center.  Coordinates go
    from -0.5 to 0.5.
    """
    if exp_type == "NRS_IFU":
        #
        # Currently assume IFU sources are centered
        return (0.0, 0.0)
    elif exp_type == "NRS_MSASPEC":
        #
        # MSA centering specified in the MiltiSlit model
        # "input" treated as a slit object
        try:
            xcenter = input.source_xcent
            ycenter = input.source_ycent
        except AttributeError:
            log.warning("Unable to get source center from model")
            log.warning("Using 0.0, 0.0")
            xcenter = 0.0
            ycenter = 0.0
        return (xcenter, ycenter)
    else:
        log.warning("No method to get centering for exp_type %s" % exp_type)
        log.warning("Using (0.0, 0.0)")
        return (0.0, 0.0)

def getApertureFromModel(input_model, nshutters):
    """Figure out the correct aperture based on the number of shutters in the slit
    """
    for aperture in input_model.apertures:
        if aperture.shutters == nshutters: return aperture
    #
    # If nothing matches, return None
    return None

def calculate_pathloss_vector(pathloss_model_type, xcenter, ycenter):
    """
    Calculate the pathloss vector from the pathloss model using the
    coordinates of the center of the target to interpolate the
    pathloss value as a function of wavelength at that location

    Parameters:
    -----------

    pathloss_model_type: datamodel.PathlossModel.pointsource
                      or datamodel.PathlossModel.uniformsource

    The input pathloss model attribute

    xcenter: Float

    The x-center of the target (-0.5 to 0.5)

    ycenter: Float

    The y-center of the target (-0.5 to 0.5)

    """

    refdata = pathloss_model_type.data
    pathwcs = pathloss_model_type.wcs
    wavesize = refdata.shape[0]
    wavelength = np.zeros(wavesize, dtype=np.float32)
    #
    # uniformsource.data is 1-d, we just return it, along with
    # a vector of wavelengths calculated using the WCS
    if len(refdata.shape) == 1:
        crpix1 = pathwcs.crpix1
        crval1 = pathwcs.crval1
        cdelt1 = pathwcs.cdelt1    
        for i in np.arange(wavesize):
            wavelength[i] = crval1 +(float(i) - crpix1)*cdelt1
        return wavelength, refdata
    #
    # pointsource.data is 3-d, so we have to extract a wavelength vector
    # at the specified location.  We do this using bilinear interpolation
    else:
        crpix3 = pathwcs.crpix3
        crval3 = pathwcs.crval3
        cdelt3 = pathwcs.cdelt3    
        for i in np.arange(wavesize):
            wavelength[i] = crval3 +(float(i) - crpix3)*cdelt3
        # Calculate python index of object center
        crpix1 = pathwcs.crpix1
        crval1 = pathwcs.crval1
        cdelt1 = pathwcs.cdelt1
        crpix2 = pathwcs.crpix2
        crval2 = pathwcs.crval2
        cdelt2 = pathwcs.cdelt2
        object_colindex = crpix1 + (xcenter - crval1) / cdelt1 - 1
        object_rowindex = crpix2 + (ycenter - crval2) / cdelt2 - 1
        #
        # Do bilinear interpolation to get the array of path loss vs wavelength
        dx1 = object_colindex - int(object_colindex)
        dx2 = 1.0 - dx1
        dy1 = object_rowindex - int(object_rowindex)
        dy2 = 1.0 - dy1
        a11 = dx1*dy1
        a12 = dx1*dy2
        a21 = dx2*dy1
        a22 = dx2*dy2
        j, i = int(object_colindex), int(object_rowindex)
        pathloss_vector = a22*refdata[:, i, j] + a12*refdata[:, i+1, j] + \
            a21*refdata[:, i, j+1] + a11*refdata[:, i+1, j+1]
        return wavelength, pathloss_vector

def do_correction(input_model, pathloss_model):
    """
    Short Summary
    -------------
    Execute all tasks for Path Loss Correction

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    pathloss_model: pathloss model object
        pathloss correction data

    Returns
    -------
    output_model: data model object
        Science data with pathloss extensions added

    """
    exp_type = input_model.meta.exposure.type
    if exp_type == 'NRS_MSASPEC':
        slit_number = 0
        # For each slit
        for slit in input_model.slits:
            slit_number = slit_number + 1
            size = slit.data.size
            # That has data.size > 0
            if size > 0:
                # Get centering
                xcenter, ycenter = getCenter(exp_type, slit)
                # Calculate the 1-d wavelength and pathloss vectors
                # for the source position
                # Get the aperture from the reference file that matches the slit
                aperture = getApertureFromModel(pathloss_model, slit.nshutters)
                if aperture is not None:
                    wavelength_pointsource, pathloss_pointsource_vector = \
                        calculate_pathloss_vector(aperture.pointsource_data,
                                                  xcenter, ycenter)
                    wavelength_uniformsource, pathloss_uniform_vector = \
                        calculate_pathloss_vector(aperture.uniform_data,
                                                  xcenter, ycenter)
                    #
                    # Wavelengths in the reference file are in meters, need them to be
                    # in microns
                    wavelength_pointsource *= 1.0e6
                    wavelength_uniformsource *= 1.0e6
                else:
                    print("Cannot find matching pathloss model for slit with size %d" % slit.nshutters)
                    continue
                nrows, ncols = slit.data.shape
                # Get wavelengths of each end
                xstart = slit.xstart
                xstop = xstart + ncols
                ystart = slit.ystart
                ystop = ystart + nrows
                # For each pixel
                y, x = np.mgrid[ystart:ystop, xstart:xstop]
                ra, dec, wave_array = slit.meta.wcs(x, y)
                slit.pl_point = calculate_pathloss(wave_array,
                                                   wavelength_pointsource,
                                                   pathloss_pointsource_vector)
                slit.pl_uniform = calculate_pathloss(wave_array,
                                                     wavelength_uniformsource,
                                                     pathloss_uniform_vector)
    elif exp_type == 'NRS_IFU':
        #
        # Make empty pathloss arrays for point and uniform
        nrows, ncols = input_model.data.shape
        input_model.pl_point = np.ones((nrows, ncolumns), dtype=np.float32)
        input_model.pl_uniform = np.ones((nrows, ncolumns), dtype=np.float32)
        #
        # Get the list of WCS objects for each IFU slice
        wcslist = nirspec.nrs_ifu_wcs(input_model)
        for ifuslice in wcslist:
            # Get centering
            xcenter, ycenter = getCenter(exp_type, ifuslice)
            # Calculate the 1-d wavelength and pathloss vectors
            # for the source position
            wavelength_pointsource, pathloss_pointsource_vector = \
                calculate_pathloss_vector(pathloss_model.pointsource,
                                          xcenter, ycenter)
            wavelength_uniformsource, pathloss_uniform_vector = \
                calculate_pathloss_vector(pathloss_model.uniformsource,
                                          xcenter, ycenter)
            xmin = ifuslice.domain[0]['lower']
            xmax = ifuslice.domain[0]['upper']
            ymin = ifuslice.domain[1]['lower']
            ymax = ifuslice.domain[1]['upper']
            y, x = np,mgrid[ymin:ymax, xmin:xmax]
            ra, dec, wave_array = slit.meta.wcs(x, y)
            input_model.pl_point[ymin:ymax, xmin:xmax] = calculate_pathloss(wave_array,
                                                                            wavelength_pointsource,
                                                                            pathloss_pointsource_vector)
            input_model.pl_uniform[ymin:ymax, xmin:xmax] = calculate_pathloss(wave_array,
                                                                              wavelength_uniformsource,
                                                                              pathloss_uniform_vector)
            
    return input_model.copy()

def calculate_pathloss(wave_array, wavelength, pathloss):
    """
    Short Summary
    -------------
    
    Calculate pathloss given wavelength of each pixel in a 2-d array, and
    1-d arrays of wavelength and pathloss

    Parameters
    ----------

    wave_array: numpy ndarray
        Input array of wavelength values at each pixel

    wavelength: numpy ndarray
        1-d array of wavelength values

    pathloss: numpy ndarray
        1-d array of corresponding pathloss values

    Returns:
    --------

    pathloss_array: 2-d array of pathloss values

    """

    #
    # Make sure wavelength and pathloss arrays have the same length
    if wavelength.shape[0] != pathloss.shape[0]:
        log.warning("Wavelength and pathloss arrays have different dimensions")
        log.info("wavelength shape = %s, pathloss shape = %s" % (wavelength.shape[0],
                                                                 pathloss.shape[0]))
        return None

    nrows, ncols = wave_array.shape

    #
    # Make an array of ones for the pathloss
    pathloss_array = np.ones((nrows, ncols), dtype=np.float32)
    
    for i in range(nrows):
        for j in range(ncols):
            array_value = wave_array[i, j]
            pathloss_value = interpolated_lookup(array_value, wavelength,
                                                 pathloss)
            pathloss_array[i, j] = pathloss_value

    return pathloss_array

def interpolated_lookup(value, array_in, array_out):
    """
    Short Summary
    -------------
    
    Given two 1-d arrays of corresponding (x, y) values, use the input
    value to find the location in the first array using linear
    interpolation, then linearly interpolate at that location in the
    second array

    Parameters
    ----------

    value: float
        Input value the be used to calculate index in array_in

    array_in: numpy ndarray
        1-d array of input values

    array_out: numpy ndarray
        1-d array of corresponding output values

    Returns:
    --------

    output: float
        Value interpolated at the corresponding location of array_out

    """

    subtracted = array_in - value
    shift_mult = subtracted[1:] * subtracted[:-1]
    index_tuple = np.where(shift_mult[1:]*shift_mult[:-1] < 0.0)
    if len(index_tuple[0]) > 0:
        index = index_tuple[0][0]
        remainder = value - array_in[index]
        partial = remainder/(array_in[index+1] - array_in[index])
        returned = array_out[index] + \
                   partial*(array_out[index+1] - array_out[index])
        return returned
    else:
        if subtracted[0] >= 0:
            return array_out[0]
        else:
            return array_out[-1]

def get_extension(hdulist, filter):
    """Return the first extension that meets the criteria specified
    in the dictionary 'filter'

    Parameters:
    -----------
    hdulist: fits HDUList
        The fits HDUList to be searched.

    filter: dict
        A dictionary specifying the conditions to be met by an extension.  The
        conditions are a (key, value) pair with 'key' being the string containing
        the name of a FITS keyword, and 'value' being the required keyword value.
        For example, filter could be {'EXTNAME': 'SCI', 'APERTURE': 'A1600'}, and
        the routine will return the index of the first extension that matches the
        filter.

    Returns
    -------
    integer
        The index of the first extension that matches the filter.
    """
    matched_index = None
    for index, hdu in enumerate(hdulist):
        for key in filter.keys():
            if hdu.header.get(key) == filter[key]:
                matched_index = index
            else:
                matched_index = None
                break
        if matched_index is not None:
            return matched_index
    return None
