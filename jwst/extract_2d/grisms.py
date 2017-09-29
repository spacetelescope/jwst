#
#  Module for 2d extraction of grism spectra
#
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging
import copy
import numpy as np
from astropy.modeling.models import Shift, Const1D, Mapping, Identity

from .. import datamodels
from ..assign_wcs import util


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract_grism_objects(input_model, grism_objects=[], reffile=""):
    """
    Extract 2d boxes around each objects spectra per order.

    Parameters
    ----------
    input_model : jwst.datamodels.ImageModel

    grism_objects : list[GrismObject]
        A list of GrismObjects

    reffile : str
        The name of the wavelengthrange reference file


    Returns
    -------
    output_model : jwst.datamodels.MultiSlitModel()


    Notes
    -----
    This method supports NRC_GRISM and NIS_WFSS only

    GrismObject is a named tuple which contains distilled
    information about each catalog object. It can be created
    by calling jwst.assign_wcs.util.create_grism_bbox() which
    will return a list of GrismObjects that countain the bounding
    boxes needed for extraction.

    """
    if not grism_objects:
        # get the wavelengthrange reference file from the input_model
        if not reffile:
            raise ValueError("Expected name of wavelengthrange reference file")
        else:
            grism_objects = util.create_grism_bbox(input_model, {"wavelengthrange": reffile})


    output_model = datamodels.MultiSlitModel()
    output_model.update(input_model)

    # One WCS model can be used to govern all the extractions
    # and in fact the model transforms rely on the full frame
    # coordinates of the input pixel location. So the WCS
    # attached to the extraction is just a copy of the 
    # input_model WCS with a shift transform to the corner
    # of the subarray. They also depend on the source object
    # center, this information will be saved to the meta of
    # the output model as source_[x/y]pos
    inwcs = input_model.meta.wcs 
    subwcs = copy.deepcopy(inwcs)

    # For easy reference here, GrismObjects has:
    #
    # xcenter,ycenter: in direct image pixels
    # order_bounding in grism_detector pixels
    # icrs_centroid: SkyCoord of object center
    # sky_bbox_ :lower and upper bounding box in SkyCoord
    # sid: catalog ID of the object

    for obj in grism_objects:
        for order in obj.order_bounding.keys():

            # Add the shift to the lower corner to each subarray WCS object
            # The shift should just be the lower bounding box corner
            # also replace the object center location inputs to the GrismDispersion
            # model with the known object center information (in pixels of direct image)
            # This is changes the user input to the model from (x,y,x0,y0,order) -> (x,y,order)
            #
            # The bounding boxes here are also limited to the size of the detector
            bb = obj.order_bounding[order]
            xmin, xmax = int(round(max(bb[0][0], 0))), int(round(min(bb[0][1], input_model.meta.subarray.xsize)))  # limit the boxes to the detector
            ymin, ymax = int(round(max(bb[1][0], 0))), int(round(min(bb[1][1], input_model.meta.subarray.ysize))) # limit the boxes to the detector

            tr = inwcs.get_transform('grism_detector', 'detector')
            tr = Identity(3) | (Mapping((0, 1, 0, 1, 2)) | Shift(xmin) & Shift(ymin) &
                                                           Const1D(obj.xcenter) & Const1D(obj.ycenter) &
                                                           Identity(1)) | tr
            subwcs.set_transform('grism_detector', 'detector', tr)

            log.info("Subarray extracted for obj: {} order: {}:".format(obj.sid, order))
            log.info("Subarray extents are: ({}, {}), ({}, {})".format(xmin,ymin,xmax,ymax))

            ext_data = input_model.data[ymin : ymax + 1, xmin : xmax+ 1].copy()
            ext_err = input_model.err[ymin : ymax + 1, xmin : xmax + 1].copy()
            ext_dq = input_model.dq[ymin : ymax + 1, xmin : xmax + 1].copy()

            
            new_model = datamodels.ImageModel(data=ext_data, err=ext_err, dq=ext_dq)
            new_model.meta.wcs = subwcs
            output_model.slits.append(new_model)

            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            # nslit = obj.sid - 1  # catalog id starts at zero
            output_model.slits[-1].name = str(obj.sid)
            output_model.slits[-1].xstart = xmin + 1
            output_model.slits[-1].xsize = (xmax - xmin) + 1
            output_model.slits[-1].ystart = ymin + 1
            output_model.slits[-1].ysize = (ymax - ymin) + 1
        
    
    # memory reduction for pipeline chain        
    del subwcs
    return output_model


# move to extract 2d for grism exposures
# add the wavelength range using the far left and far right orders?
# fselect1 = wrange_selector[0].index(input_model.meta.instrument.filter)
# fselect2 = wrange_selector[-1].index(input_model.meta.instrument.filter)
# lower_lam = wrange[0][fselect1][0]
# upper_lam = wrange[-1][fselect2][1]
# input_model.meta.wcsinfo.waverange_start = lower_lam
# input_model.meta.wcsinfo.waverange_end = upper_lam

def compute_dispersion(wcs):
    """
    Compute the pixel dispersion.
    Make a model for the pixel dispersion from the grismconf specs

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.

    Returns
    -------
    dispersion : ndarray-like
        The pixel dispersion [in m].

    """
    raise NotImplementedError

