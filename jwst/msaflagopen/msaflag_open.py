from __future__ import division

#
#  Module for flagging the DQ array of pixels affected by failed
#  open MSA shutters in nirspec science data sets
#

import math
import numpy as np
import logging
from jwst_lib import models
from jwst_pipeline.assign_wcs import nirspec

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

BADFLAG = models.dqflags.pixel['MSA_FAILED_OPEN']

def do_correction(input_model, failedopen_model):
    """
    Short Summary
    -------------
    Apply DQ flag to pixels affected by failed open MSA shutters

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    failedopen_model: failedopen model object
        failed open reference data

    Returns
    -------
    output_model: data model object
        science data with DQ array modified

    """
       
    # Flag the stuck open shutters
    output_model = flag(input_model, failedopen)

    output_model.meta.cal_step.msaflagopen = 'COMPLETE'

    return output_model


def flag(input, failedopen):
    """
    Takes the list of failed open shutters from thje failedopen reference file
    and calculates the pixels affected using the WCS model.
    The affected pixels in the science data have their DQ flags combined with
    that for the MSA_FAILED_OPEN standard flag.  All other science data
    arrays are unchanged.

    Parameters
    ----------
    input: data model object
        the input science data

    failedopen: failedopen structure
        the failedopen reference data

    Returns
    -------
    output: data model object
        science data with DQ flags of affected modified

    """

    # Create output as a copy of the input science data model
    output = input.copy()
    filter = input.meta.instrument.filter
    grating = input.meta.instrument.grating
    filter_grating = filter + '_' + grating
    waveref = input.meta.ref_file.wavelengthrange.name
    wavetree = AsdfFile.open(waveref).tree
    wmin, wmax = wavetree['filter_grating'][filter_grating]['range']
    order = wavetree['filter_grating'][filter_grating]['order']
    slits = []
    for failed in failedopen:
        quadrant = failed['Q']
        x = failed['x']
        y = failed['y']
        slits.append((quadrant, x, y))
    slitmodel = nirspec.slit_to_detector(input, slits)
    for slit in slitmodel:
        #
        # Calculate length of slit
        # Length of slit at minimum wavelength
        x1, y1 = slit(0.0, -0.5, wmin)
        x2, y2 = slit(0.0, 0.5, wmin)
        minlength = math.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        # Length at maximum wavelength
        x1, y1 = slit(0.0, -0.5, wmax)
        x2, y2 = slit(0.0, 0.5, wmax)
        maxlength = math.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        # Average these to get the slit length
        slitlength = 0.5*(maxlength + minlength)
        setdq(output, xmin, ymin, xmax, ymax,
              slitlength)
    return output

def setdq(datamodel, xmin, ymin, xmax, ymax,
          slitlength):
    #
    # Set the dq data to dqflags['MSAOPEN'] when it's inside the
    # slit trace
    nrows, ncols = datamodel.dq.shape
    #
    # Calculate slope of trace
    xmin, ymin = slit(0.0, 0.0, wmin)
    xmax, ymax = slit(0.0, 0.0, wmax)
    slope = (ymax - ymin) / (xmax - xmin)
    #
    # Calculate vertical extent of slit using simple geometry
    costheta = sqrt(1.0/(1 + slope**2))
    halfheight = slitlength / (2.0*costheta)
    #
    # Round up
    halfheight = int(halfheight) + 1
    if xmin < ncols and xmax > 0:
        ixmin = max(0, int(xmin))
        ixmax = min(ncols, int(xmax))
        for ix in range(ixmin, ixmax+1):
            y = ymin + slope*(ix - xmin)
            #
            # Nearest integer
            iy = int(y + 0.5)
            dataslice = datamodel[dq][iy-halfheight:iy+halfheight+1, ix]
            dataslice = np.bitwise_or(dataslice, BADFLAG)
    return
