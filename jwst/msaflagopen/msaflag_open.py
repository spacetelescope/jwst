"""Flag pixels affected by open MSA shutters in NIRSpec exposures."""

import json
import logging
import warnings
from pathlib import Path

import numpy as np
from gwcs.wcs import WCS
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms.models import Slit

from jwst.assign_wcs.nirspec import (
    generate_compound_bbox,
    nrs_wcs_set_input,
    slitlets_wcs,
    log as nirspec_log,
)
from jwst.lib.basic_utils import LoggingContext

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


FAILEDOPENFLAG = datamodels.dqflags.pixel["MSA_FAILED_OPEN"]
SHUTTERS_PER_ROW = 365

# States in the msaoper file that are flagged when set to 'open'
FLAGGABLE_STATES = ["Internal state", "TA state", "state"]


def do_correction(input_model, shutter_refname, wcs_refnames):
    """
    Apply DQ flag to pixels affected by failed open MSA shutters.

    Parameters
    ----------
    input_model : DataModel
        Science data to be corrected.
    shutter_refname : str
        Name of MSAOPER reference file.
    wcs_refnames : dict
        Dictionary of wcs reference file names.

    Returns
    -------
    output_model : DataModel
        Science data with DQ array modified.
    """
    # Create a list of failed open slitlets from the msaoper reference file
    failed_slitlets = create_slitlets(shutter_refname)
    log.info("%d failed open shutters", len(failed_slitlets))

    # Flag the stuck open shutters
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="Invalid interval")
        output_model = flag(input_model, failed_slitlets, wcs_refnames)
    output_model.meta.cal_step.msa_flagging = "COMPLETE"

    return output_model


def flag(input_datamodel, failed_slitlets, wcs_refnames):
    """
    Flag slitlet regions for failed open shutters.

    Takes the list of failed open shutters from the failedopen reference file
    and calculates the pixels affected using the WCS model.
    The affected pixels in the science data have their DQ flags combined with
    that for the MSA_FAILED_OPEN standard flag.  All other science data
    arrays are unchanged.

    The input datamodel is modified in place.

    Parameters
    ----------
    input_datamodel : DataModel
        Input science data.
    failed_slitlets : list of Slit
        Failed open slitlets.
    wcs_refnames : dict
        Reference file names used to calculate the WCS. Keys are reference
        file types; values are file paths.

    Returns
    -------
    DataModel
        Science data with DQ flags modified.
    """
    # Use the machinery in assign_wcs to create a WCS object for the bad shutters
    with LoggingContext(nirspec_log, level=logging.WARNING):
        pipeline = slitlets_wcs(input_datamodel, wcs_refnames, failed_slitlets)
    wcs = WCS(pipeline)

    # Create output as a copy of the input science data model so we can overwrite
    # the wcs with the wcs for the failed open shutters
    # Have to make sure the EXP_TYPE is NRS_MSASPEC so that nrs_wcs_set_input works,
    # We need to use the slit WCS for this even if the EXP_TYPE is NRS_IFU because we
    # are calculating where stuck open slits affect the data
    temporary_copy = input_datamodel.copy()
    temporary_copy.meta.wcs = wcs
    temporary_copy.meta.exposure.type = "NRS_MSASPEC"
    temporary_copy.meta.wcs.bounding_box = generate_compound_bbox(temporary_copy, failed_slitlets)

    dq_array = input_datamodel.dq
    for slitlet in failed_slitlets:
        # Pick the WCS for this slitlet from the WCS of the exposure
        thiswcs = nrs_wcs_set_input(temporary_copy, slitlet.name)

        # Convert the bounding box for this slitlet to a set of indices to use as a slice
        xmin, xmax, ymin, ymax = boundingbox_to_indices(temporary_copy, thiswcs.bounding_box)

        # Make a grid of points within the slice
        y_indices, x_indices = np.mgrid[ymin:ymax, xmin:xmax]

        # Calculate the arrays of coordinates for each pixel in the slice
        coordinate_array = thiswcs(x_indices, y_indices)

        # The coordinate_array is a tuple of arrays, one for each output coordinate
        # In this case there should be 3 arrays, one each for RA, Dec and Wavelength
        # For pixels outside the slitlet, the arrays have NaN

        # Make a subarray from these coordinate arrays by setting pixels that aren't
        # NaN to FAILEDOPENFLAG, the rest to 0
        dq_subarray = wcs_to_dq(coordinate_array, FAILEDOPENFLAG)

        # Bitwise-or this subarray with the slice in the original exposure's DQ array
        dq_array[..., ymin:ymax, xmin:xmax] |= dq_subarray

    # Set the dq array of the input datamodel to the corrected dq array
    input_datamodel.dq = dq_array
    return input_datamodel


def boundingbox_to_indices(data_model, bounding_box):
    """
    Translate a bounding box to image indices.

    Takes a bounding_box (tuple of tuples: ((x1, x2), (y1, y2)) and
    a datamodel and calculates the range of indices in the X and Y dimensions
    of the overlap between the bounding box and the datamodel's data array.

    Parameters
    ----------
    data_model : DataModel
        The input science datamodel.
    bounding_box : tuple of tuple
        Bounding box returned from wcs object.

    Returns
    -------
    xmin, xmax, ymin, ymax : int
        Range of indices of overlap between science data array and bounding box.
    """
    nrows, ncols = data_model.data.shape[-2:]
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


def wcs_to_dq(wcs_array, flag):
    """
    Create a DQ subarray corresponding to a failed open slitlet.

    The created array has the value `flag` wherever the WCS coordinates
    are valid (non-NaN) and 0 otherwise.

    Parameters
    ----------
    wcs_array : tuple of ndarray
        Image coordinates for the failed open region.
    flag : int
        DQ flag to set.

    Returns
    -------
    dq : ndarray of int
        Output DQ array.
    """
    dq = np.zeros((wcs_array[0].shape), dtype=np.uint32)
    non_nan = ~np.isnan(wcs_array[0])
    dq[non_nan] = flag
    return dq


def get_failed_open_shutters(shutter_refname):
    """
    Get the failed open shutters from a reference file.

    Parameters
    ----------
    shutter_refname : str
        File name for the MSAOPER reference file.

    Returns
    -------
    failedopen : list
        A list of shutters which satisfy the condition that at
        least one of the states in FLAGGABLE_STATES is set to 'open'.
    """
    # Read the bad shutter reference file data model
    with Path(shutter_refname).open("r") as f1:
        shutters = json.load(f1)

    failedopen = []
    for shutter in shutters["msaoper"]:
        for state in FLAGGABLE_STATES:
            if shutter[state] == "open":
                failedopen.append(shutter)
                break
    return failedopen


def create_slitlets(shutter_refname):
    """
    Create slitlets for each failed open shutter.

    For the created slit objects, "shutter_id" is an integer that uniquely defines
    the shutter in the quadrant, calculated from the x and y center of the shutter.

    In the Slit tuple, the only values that matter are "name" (must be unique),
    xcen, ycen, quadrant (from msaoper file), ymin, ymax (should be -0.5, 0.5), and
    shutter state (should be 'x', for one open shutter).  Default values are assigned
    for all other values.

    Returns
    -------
    slitlets : list of Slit
        A list of slitlets. Each slitlet is a named tuple with elements
        ("name", "shutter_id", "dither_position", "xcen", "ycen", "ymin", "ymax",
        "quadrant", "source_id", "shutter_state", "source_name", "source_alias",
        "stellarity", "source_xpos", "source_ypos", "source_ra", "source_dec").
    """
    failedopenlist = get_failed_open_shutters(shutter_refname)

    slitlets = []
    counter = 0
    for shutter in failedopenlist:
        counter = counter + 1
        x = shutter["x"]
        y = shutter["y"]
        shutter_id = x + (y - 1) * SHUTTERS_PER_ROW
        slitlets.append(
            Slit(counter, shutter_id, 0, x, y, -0.5, 0.5, shutter["Q"], 0, "x", slit_id=counter)
        )
    return slitlets
