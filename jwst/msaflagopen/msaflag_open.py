#
#  Module for flagging the DQ array of pixels affected by failed
#  open MSA shutters in nirspec science data sets
#
import json
import numpy as np
import logging
from .. import datamodels
from ..assign_wcs.nirspec import slitlets_wcs, nrs_wcs_set_input
from ..transforms.models import Slit
from gwcs.wcs import WCS

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


FAILEDOPENFLAG = datamodels.dqflags.pixel['MSA_FAILED_OPEN']
SHUTTERS_PER_ROW = 365
#
# States in the msaoper file that are to flagged when set to 'open'
FLAGGABLE_STATES = ['Internal state', 'TA state', 'state']


def do_correction(input_model, shutter_refname, wcs_refnames):
    """
    Short Summary
    -------------
    Apply DQ flag to pixels affected by failed open MSA shutters

    Parameters
    ----------
    input_model: datamodel object
        science data to be corrected

    shutter_refname: string
        Name of MSAOPER reference file

    wcs_refnames: dict
        Dictionary of wcs reference file names

    Returns
    -------
    output_model: data model object
        science data with DQ array modified

    """
    #
    # Create a list of failed open slitlets from the msaoper reference file
    failed_slitlets = create_slitlets(input_model, shutter_refname)

    # Flag the stuck open hutters
    output_model = flag(input_model, failed_slitlets, wcs_refnames)

    output_model.meta.cal_step.msa_flagging = 'COMPLETE'

    return output_model


def flag(input_datamodel, failed_slitlets, wcs_refnames):
    """
    Takes the list of failed open shutters from the failedopen reference file
    and calculates the pixels affected using the WCS model.
    The affected pixels in the science data have their DQ flags combined with
    that for the MSA_FAILED_OPEN standard flag.  All other science data
    arrays are unchanged.

    Parameters
    ----------
    input_datamodel: data model object
        the input science data

    failed_slitlets: list
        List of failed open slitlets

    wcs_refnames: dict
        dictionary of reference file names used to calculate the WCS

    Returns
    -------
    input_datamodel: data model object
        science data with DQ flags of affected modified

    """
    #
    # Use the machinery in assign_wcs to create a WCS object for the bad shutters
    pipeline = slitlets_wcs(input_datamodel, wcs_refnames, failed_slitlets)
    wcs = WCS(pipeline)
    # Create output as a copy of the input science data model so we can overwrite
    # the wcs with the wcs for the failed open shutters
    # Have to make sure the EXP_TYPE is NRS_MSASPEC so that nrs_wcs_set_input works,
    # We need to use the slit WCS for this even if the EXP_TYPE is NRS_IFU because we
    # are calculating where stuck open slits affect the data
    temporary_copy = input_datamodel.copy()
    temporary_copy.meta.wcs = wcs
    temporary_copy.meta.exposure.type = 'NRS_MSASPEC'
    dq_array = input_datamodel.dq
    for slitlet in failed_slitlets:
        #
        # Pick the WCS for this slitlet from the WCS of the exposure
        thiswcs = nrs_wcs_set_input(temporary_copy, slitlet.name)
        #
        # Convert the bounding box for this slitlet to a set of indices to use as a slice
        xmin, xmax, ymin, ymax = boundingbox_to_indices(temporary_copy,
                                                        thiswcs.bounding_box)
        #
        # Make a grid of points within the slice
        y_indices, x_indices = np.mgrid[ymin:ymax, xmin:xmax]
        #
        # Calculate the arrays of coordinates for each pixel in the slice
        coordinate_array = thiswcs(x_indices, y_indices)
        #
        # The coordinate_array is a tuple of arrays, one for each output coordinate
        # In this case there should be 3 arrays, one each for RA, Dec and Wavelength
        # For pixels outside the slitlet, the arrays have NaN
        #
        # Make a subarray from these coordinate arrays by setting pixels that aren't
        # NaN to FAILDOPENFLAG, the rest to 0
        dq_subarray = wcs_to_dq(coordinate_array, FAILEDOPENFLAG)
        #
        # bitwise-or this subarray with the slice in the original exposure's DQ array
        dq_array = or_subarray_with_array(dq_array, dq_subarray, xmin, xmax, ymin, ymax)
    #
    # Set the dq array of the input datamodel to the corrected dq array
    input_datamodel.dq = dq_array
    return input_datamodel


def boundingbox_to_indices(data_model, bounding_box):
    """
    Takes a bounding_box (tuple of tuples: ((x1, x2), (y1, y2)) and
    a datamodel and calculates the range of indices in the X and Y dimensions
    of the overlap between the bounding box and the datamodel's data array.

    Parameters
    ----------
    data_model: data model object
        the input science datamodel

    bounding_box: Tuple of tuples
        Bounding box returned from wcs object

    Returns
    -------
    xmin, xmax, ymin, ymax: integers
        Range of indices of overlap between science data array and bounding box

    """
    nrows, ncols = data_model.data.shape
    ((x1, x2), (y1, y2)) = bounding_box
    xmin = int(min(x1, x2))
    xmin = max(xmin, 0)
    xmax = int(max(x1, x2)) + 1
    xmax = min(xmax, ncols)
    ymin = int(min(y1, y2))
    ymin = max(ymin, 0)
    ymax = int(max(y1, y2)) + 1
    ymax = min(ymax, nrows)
    return xmin, xmax, ymin, ymax


def wcs_to_dq(wcs_array, FLAG):
    #
    # Convert one of the arrays in the wcs_array tuple to a DQ array,
    # which has zeros everywhere the wcs array has the value NaN, and
    # has the value FLAG everywhere the wcs array is not NaN
    dq = np.zeros((wcs_array[0].shape), dtype=np.uint32)
    non_nan = np.where(~np.isnan(wcs_array[0]))
    dq[non_nan] = FLAG
    return dq


def get_failed_open_shutters(shutter_refname):
    """
    Return a list of shutters which satisfy the condition that at
    least one of the states in FLAGGABLE_STATES is set to 'open'
    """
    # Open the bad shutter reference file data model
    f1 = open(shutter_refname, 'r')
    shutters = json.load(f1)
    f1.close()

    failedopen = []
    for shutter in shutters['msaoper']:
        for state in FLAGGABLE_STATES:
            if shutter[state] == 'open':
                failedopen.append(shutter)
                break
    return failedopen


def create_slitlets(input_model, shutter_refname):
    """
    A slitlet looks like this:

    slitlets : list
        A list of slitlets. Each slitlet is a named tuple with
        ("name", "shutter_id", "xcen", "ycen", "ymin", "ymax", "quadrant",
        "source_id", "shutter_state")

    A slit is:

    Slit = namedtuple('Slit', ["name", "shutter_id", "dither_position", "xcen", "ycen",
                           "ymin", "ymax", "quadrant", "source_id", "shutter_state",
                           "source_name", "source_alias", "stellarity",
                           "source_xpos", "source_ypos"])
    "shutter_id" is an integer that uniquely defines the shutter in the quadrant, it is calculated
    from the x and y using the function
    Slit.__new__.__defaults__= ("", 0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, "", "", "", 0.0, 0.0, 0.0)

    The only ones that matter are "name" (must be unique), xcen, ycen, quadrant (from msaoper
    file), ymin, ymax (should be -0.5, 0.5), nshutters (should be 1)
    """

    failedopenlist = get_failed_open_shutters(shutter_refname)

    slitlets = []
    counter = 0
    for shutter in failedopenlist:
        counter = counter + 1
        x = shutter['x']
        y = shutter['y']
        shutter_id = id_from_xy(x, y)
        slitlets.append(Slit(str(counter), shutter_id, 0, x, y, -0.5, 0.5,
                             shutter['Q'], 0, 1, "", "", 0.0, 0.0, 0.0)
                        )
    return slitlets


def id_from_xy(x, y):
    """
    Calculate the shutter_id from the x and y of a shutter
    """

    shutter_id = x + (y - 1) * SHUTTERS_PER_ROW
    return shutter_id


def or_subarray_with_array(dq_array, dq_subarray, xmin, xmax, ymin, ymax):
    """
    Bitwise-or the slice of the dq array with the section corresponding to the
    failed open shutter
    """
    dq_array[ymin:ymax, xmin:xmax] = np.bitwise_or(dq_array[ymin:ymax, xmin:xmax], dq_subarray)
    return dq_array
