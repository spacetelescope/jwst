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
    elif exp_type in ["NRS_MSASPEC", "NRS_FIXEDSLIT", "NRS_BRIGHTOBJ"]:
        #
        # MSA centering specified in the MiltiSlit model
        # "input" treated as a slit object
        try:
            xcenter = input.source_xpos
            ycenter = input.source_ypos
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

def getApertureFromModel(input_model, match):
    """Figure out the correct aperture based on the value of the 'match'
    parameter.  For MSA, match is the number of shutters, for fixed slit,
    match is the name of the slit
    """
    if input_model.meta.exposure.type == 'NRS_MSASPEC':
        for aperture in input_model.apertures:
            if aperture.shutters == match: return aperture
    elif input_model.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
        for aperture in input_model.apertures:
            log.debug(aperture.name)
            if aperture.name == match: return aperture
    else:
        log.warning('Unable to get aperture from model type {0}'.format(input_model.meta.exposure.type))
    #
    # If nothing matches, return None
    return None

def calculate_pathloss_vector(pathloss_refdata, pathloss_wcs, xcenter, ycenter):
    """
    Calculate the pathloss vectors from the pathloss model using the
    coordinates of the center of the target to interpolate the
    pathloss value as a function of wavelength at that location

    Parameters:
    -----------

    pathloss_refdata:     numpy ndarray

    The input pathloss data array

    pathloss_wcs:      wcs attribute from model

    xcenter: Float

    The x-center of the target (-0.5 to 0.5)

    ycenter: Float

    The y-center of the target (-0.5 to 0.5)

    """
    
    wavesize = pathloss_refdata.shape[0]
    wavelength = np.zeros(wavesize, dtype=np.float32)
    #
    # uniformsource.data is 1-d, we just return it, along with
    # a vector of wavelengths calculated using the WCS
    if len(pathloss_refdata.shape) == 1:
        crpix1 = pathloss_wcs.crpix1
        crval1 = pathloss_wcs.crval1
        cdelt1 = pathloss_wcs.cdelt1    
        for i in np.arange(wavesize):
            wavelength[i] = crval1 +(float(i) - crpix1)*cdelt1
        return wavelength, pathloss_refdata
    #
    # pointsource.data is 3-d, so we have to extract a wavelength vector
    # at the specified location.  We do this using bilinear interpolation
    else:
        crpix3 = pathloss_wcs.crpix3
        crval3 = pathloss_wcs.crval3
        cdelt3 = pathloss_wcs.cdelt3    
        for i in np.arange(wavesize):
            wavelength[i] = crval3 +(float(i) - crpix3)*cdelt3
        # Calculate python index of object center
        crpix1 = pathloss_wcs.crpix1
        crval1 = pathloss_wcs.crval1
        cdelt1 = pathloss_wcs.cdelt1
        crpix2 = pathloss_wcs.crpix2
        crval2 = pathloss_wcs.crval2
        cdelt2 = pathloss_wcs.cdelt2
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
        pathloss_vector = a22*pathloss_refdata[:, i, j] + a12*pathloss_refdata[:, i+1, j] + \
            a21*pathloss_refdata[:, i, j+1] + a11*pathloss_refdata[:, i+1, j+1]
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
    log.info(exp_type)
    if exp_type == 'NRS_MSASPEC':
        slit_number = 0
        # For each slit
        for slit in input_model.slits:
            slit_number = slit_number + 1
            log.info('Working on slit %d' % slit_number)
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
                                                  aperture.pointsource_wcs,
                                                  xcenter, ycenter)
                    wavelength_uniformsource, pathloss_uniform_vector = \
                        calculate_pathloss_vector(aperture.uniform_data,
                                                  aperture.uniform_wcs,
                                                  xcenter, ycenter)
                    #
                    # Wavelengths in the reference file are in meters, need them to be
                    # in microns
                    wavelength_pointsource *= 1.0e6
                    wavelength_uniformsource *= 1.0e6
                    slit.pathloss_pointsource = pathloss_pointsource_vector
                    slit.wavelength_pointsource =  wavelength_pointsource
                    slit.pathloss_uniformsource = pathloss_uniform_vector
                    slit.wavelength_uniformsource = wavelength_uniformsource
                else:
                    log.warning("Cannot find matching pathloss model for slit with size %d" % slit.nshutters)
                    continue
        input_model.meta.cal_step.pathloss = 'COMPLETE'
    elif exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
        slit_number = 0
        # For each slit
        for slit in input_model.slits:
            log.info(slit.name)
            slit_number = slit_number + 1
            # Get centering
            xcenter, ycenter = getCenter(exp_type, slit)
            # Calculate the 1-d wavelength and pathloss vectors
            # for the source position
            # Get the aperture from the reference file that matches the slit
            aperture = getApertureFromModel(pathloss_model, slit.name)
            if aperture is not None:
                log.info("Using aperture {0}".format(aperture.name))
                wavelength_pointsource, pathloss_pointsource_vector = \
                    calculate_pathloss_vector(aperture.pointsource_data,
                                              aperture.pointsource_wcs,
                                              xcenter, ycenter)
                wavelength_uniformsource, pathloss_uniform_vector = \
                    calculate_pathloss_vector(aperture.uniform_data,
                                              aperture.uniform_wcs,
                                              xcenter, ycenter)
                #
                # Wavelengths in the reference file are in meters, need them to be
                # in microns
                wavelength_pointsource *= 1.0e6
                wavelength_uniformsource *= 1.0e6
                slit.pathloss_pointsource = pathloss_pointsource_vector
                slit.wavelength_pointsource =  wavelength_pointsource
                slit.pathloss_uniformsource = pathloss_uniform_vector
                slit.wavelength_uniformsource = wavelength_uniformsource
            else:
                log.warning("Cannot find matching pathloss model for aperture %s" % slit.name)
                continue
        input_model.meta.cal_step.pathloss = 'COMPLETE'
    elif exp_type == 'NRS_IFU':
        # Get centering
        xcenter, ycenter = getCenter(exp_type, None)
        # Calculate the 1-d wavelength and pathloss vectors
        # for the source position
        aperture = pathloss_model.apertures[0]
        wavelength_pointsource, pathloss_pointsource_vector = \
            calculate_pathloss_vector(aperture.pointsource_data,
                                      aperture.pointsource_wcs,
                                      xcenter, ycenter)
        wavelength_uniformsource, pathloss_uniform_vector = \
            calculate_pathloss_vector(aperture.uniform_data,
                                      aperture.uniform_wcs,
                                      xcenter, ycenter)
        # Wavelengths in the reference file are in meters, need them to be
        # in microns
        wavelength_pointsource *= 1.0e6
        wavelength_uniformsource *= 1.0e6
        input_model.wavelength_pointsource = wavelength_pointsource
        input_model.pathloss_pointsource = pathloss_pointsource_vector
        input_model.wavelength_uniformsource = wavelength_uniformsource
        input_model.pathloss_uniformsource = pathloss_uniform_vector
        input_model.meta.cal_step.pathloss = 'COMPLETE'
            
    return input_model.copy()
