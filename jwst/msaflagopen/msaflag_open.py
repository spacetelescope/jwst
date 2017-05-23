from __future__ import division

#
#  Module for flagging the DQ array of pixels affected by failed
#  open MSA shutters in nirspec science data sets
#

import math
import numpy as np
import logging
from .. import datamodels
from ..assign_wcs import nirspec
from ..transforms.models import Slit

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


BADFLAG = datamodels.dqflags.pixel['MSA_FAILED_OPEN']

FLAGGABLE_STATES = ['Internal state', 'TA state', 'state']

def do_correction(input_model, shutters):
    """
    Short Summary
    -------------
    Apply DQ flag to pixels affected by failed open MSA shutters

    Parameters
    ----------
    input_model: datamodel object
        science data to be corrected

    shutters: list
        list of shutters

    Returns
    -------
    output_model: data model object
        science data with DQ array modified

    """
    failedopenlist = get_failed_open_shutters(shutters, flaggable_states)

    # Flag the stuck open shutters
    output_model = flag(input_model, failedopenlist)

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
    #
    # Use Nadia's machinery to create a WCS object for each bad shutter
    failed_slitlets = create_slitlets(faildopen)
    pipeline = create_pipeline(input_model, failed_slitlets)
    wcs = WCS(pipeline)
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

def get_failed_open_shutters(shutters):
    """
    Return a list of shutters which satisfy the condition that at
    least one of the states in FLAGGABLE_STATES is set to 'open'
    """
    failedopen = []
    for shutter in shutters:
        for state in FLAGGABLE_STATES:
            if shutter[state] == 'open':
                failedopen.append(shutter)
                break
    return failedopen

def create_slitlets(failedopen):
    """A slitlet looks like this:
    slitlets : list
        A list of slitlets. Each slitlet is a named tuple with
        ("name", "shutter_id", "xcen", "ycen", "ymin", "ymax", "quadrant", "source_id", "nshutters")

    A slit is:
    
    Slit = namedtuple('Slit', ["name", "shutter_id", "xcen", "ycen",
                           "ymin", "ymax", "quadrant", "source_id", "nshutters",
                           "source_name", "source_alias", "catalog_id", "stellarity",
                           "source_xpos", "source_ypos"])
    Slit.__new__.__defaults__= ("", 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, "", "", "", 0.0, 0.0, 0.0)

    The only ones that matter are "name" (must be unique), xcen, ycen, quadrant (from msaoper
    file), ymin, ymax (should be -0.5, 0.5), nshutters (should be 1)
 
    """
    slitlets = []
    counter = 0
    for shutter in failedopen:
        counter = counter + 1
        slitlets.append(Slit(str(counter), 0, shutter['x'], shutter['y'], -0.5, 0.5,
                             shutter['Q'], 0, 1, "", "", "", 0.0, 0.0, 0.0)
                        )
    return slitlets
