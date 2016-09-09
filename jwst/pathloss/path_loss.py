from __future__ import division

#
#  Module for calculating pathloss correction for science data sets
#

import numpy as np
import logging
from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def getCenter(input_model):
    """
    Get the center of the target in the aperture.
    (0.0, 0.0) is the aperture center.  Coordinates go
    from -0.5 to 0.5.
    """
    return (0.0, 0.0)

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

    # Get centering
    x0, y0 = getCenter(input_model)
    refdata = pathloss_model.pointsource.data
    wavesize, ysize, xsize = refdata.shape
    pathwcs = pathloss_model.pointsource.wcs
    crpix1 = pathwcs.crpix1
    crval1 = pathwcs.crval1
    cdelt1 = pathwcs.cdelt1
    crpix2 = pathwcs.crpix2
    crval2 = pathwcs.crval2
    cdelt2 = pathwcs.cdelt2
    crpix3 = pathwcs.crpix3
    crval3 = pathwcs.crval3
    cdelt3 = pathwcs.cdelt3    
    # Calculate python index of object center
    object_colindex = crpix1 + (x0 - crval1) / cdelt1 - 1
    object_rowindex = crpix2 + (y0 - crval2) / cdelt2 - 1
    print("Model spatial extent is %d by %d pixels" % (xsize, ysize))
    print("Object center is %d, %d (0-indexed)" % (object_colindex,
                                                   object_rowindex))
    #
    # Do bilinear interpolation to get the array of path loss vs wavelength
    wavelength = np.zeros(wavesize, dtype=np.float32)
    for i in np.arange(wavesize).astype(np.float32):
        wavelength[i] = crval3 +(i-crpix3)*cdelt3
    dx1 = object_colindex - int(object_colindex)
    dx2 = 1.0 - dx1
    dy1 = object_rowindex - int(object_rowindex)
    dy2 = 1.0 - dy1
    a11 = dx1*dy1
    a12 = dx1*dy2
    a21 = dx2*dy1
    a22 = dx2*dy2
    j, i = int(object_colindex), int(object_rowindex)
    pathloss = a22*refdata[j, i, :] + a12*refdata[j+1, i, :] + \
               a21*refdata[j, i+1, :] + a11*refdata[j+1, i+1, :]
    slit_number = 0
    # For each slit
    for slit in input_model.slits:
        slit_number = slit_number + 1
        size = slit.data.size
        # That has data.size > 0
        if size > 0:
            nrows, ncols = slit.data.shape
            # Create pathloss arrays
            pathloss_point = np.zeros((nrows, ncols), dtype=np.float32)
            pathloss_uniform = np.zeros((nrows, ncols), dtype=np.float32)
            wavelength = np.zeros((nrows, ncols), dtype=np.float32)
            # Get wavelengths of each end
            xstart = slit.xstart
            xstop = xstart + slit.xsize - 1
            ystart = slit.ystart
            ystop = ystart + slit.ysize - 1
            ycenter = 0.5*(ystart + ystop)
            xmin, ymin, min_wavelength = slit.meta.wcs(xstart, ycenter)
            xmax, ymax, max_wavelength = slit.meta.wcs(xstop, ycenter)
            print("Slit %d has wavelength limits %f, %f" % (slit_number,
                                                            min_wavelength,
                                                            max_wavelength))
            print("Slit %d has parameters %d %d %d %d" % (slit_number,
                                                          xstart, xstop,
                                                          ystart, ystop))
            # For each pixel
            # Use linear interpolation as using the WCS is too slow
            x, y = np.mgrid[xstart:xstop, ystart:ystop]
            pathloss_point = calculate_pathloss(wave_array, wavelength,
                                                pathloss)

            
            
    # Get wavelength
    # Interpolate path loss
    # Set pixel in pathloss array to interpolated value
    # Write out data

    # For each SCI extension in input
    # Save some data params for easy use later

    return pathloss

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
        remainder = value - array_in[value]
        partial = remainder/(array_in[index+1] - array_in[index])
        returned = array_out[index] + \
                   partial*(array_out[index+1] - array_out[index])
