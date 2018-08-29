#
#  Module for 2d extraction of grism spectra
#

import logging
import copy
from astropy.modeling.models import Shift, Const1D, Mapping, Identity

from .. import datamodels
from ..assign_wcs import util


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract_grism_objects(input_model, grism_objects=[], reference_files={}):
    """
    Extract 2d boxes around each objects spectra per order.

    Parameters
    ----------
    input_model : jwst.datamodels.ImageModel

    grism_objects : list(GrismObject)
        A list of GrismObjects

    reffile : dict
        Needs to include the name of the wavelengthrange reference file


    Returns
    -------
    output_model : jwst.datamodels.MultiSlitModel()


    Notes
    -----
    This method supports NRC_WFSS and NIS_WFSS only

    GrismObject is a named tuple which contains distilled
    information about each catalog object. It can be created
    by calling jwst.assign_wcs.util.create_grism_bbox() which
    will return a list of GrismObjects that countain the bounding
    boxes that will be used to define the 2d extraction area.

    """
    if not grism_objects:
        # get the wavelengthrange reference file from the input_model
        if (not reference_files['wavelengthrange'] or reference_files['wavelengthrange'] == 'N/A'):
            raise ValueError("Expected name of wavelengthrange reference file")
        else:
            grism_objects = util.create_grism_bbox(input_model, reference_files)
            log.info("Grism object list created from source catalog: {0:s}"
                     .format(input_model.meta.source_catalog.filename))

    if not isinstance(grism_objects, list):
            raise ValueError("Expected input grism objects to be a list")

    log.info("Creating output model")
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
    # sky_centroid: SkyCoord of object center
    # sky_bbox_ :lower and upper bounding box in SkyCoord
    # sid: catalog ID of the object

    slits = []
    for obj in grism_objects:
        for order in obj.order_bounding.keys():

            # Add the shift to the lower corner to each subarray WCS object
            # The shift should just be the lower bounding box corner
            # also replace the object center location inputs to the GrismDispersion
            # model with the known object center and order information (in pixels of direct image)
            # This is changes the user input to the model from (x,y,x0,y0,order) -> (x,y)
            #
            # The bounding boxes here are also limited to the size of the detector
            # The check for boxes entirely off the detector is done in create_grism_bbox right now
            bb = obj.order_bounding[order]

            # limit the boxes to the detector
            ymin, ymax = (max(bb[0][0], 0), min(bb[0][1], input_model.meta.subarray.ysize))
            xmin, xmax = (max(bb[1][0], 0), min(bb[1][1], input_model.meta.subarray.xsize))

            # only the first two numbers in the Mapping are used
            # the order and source position are put directly into
            # the new wcs for the subarray for the forward transform
            tr = inwcs.get_transform('grism_detector', 'detector')
            tr = Identity(2) | Mapping((0, 1, 0, 1, 0)) | (Shift(xmin) & Shift(ymin) &
                                                           Const1D(obj.xcentroid) &
                                                           Const1D(obj.ycentroid) &
                                                           Const1D(order)) | tr

            subwcs.set_transform('grism_detector', 'detector', tr)

            # now do the same thing for the backwards transform
            # this sends wavelength along with known x0, y0, order
            # This needs to update the distortion, v2v3 transforms as well?
            # figure out how to match the inputs and outputs correctly
            # going this direction
            # tr = inwcs.get_transform('detector', 'grism_detector')
            # tr = Identity(4) | Mapping((0, 1, 2, 3)) | (Const1D(obj.xcentroid) &
            #                                             Const1D(obj.ycentroid) &
            #                                             Identity(1) &
            #                                             Const1D(order)) | tr

            # subwcs.set_transform('detector', 'grism_detector', tr)

            log.info("Subarray extracted for obj: {} order: {}:".format(obj.sid, order))
            log.info("Subarray extents are: (xmin:{}, ymin:{}), (xmax:{}, ymax:{})".format(xmin,ymin,xmax,ymax))

            ext_data = input_model.data[ymin : ymax + 1, xmin : xmax + 1].copy()
            ext_err = input_model.err[ymin : ymax + 1, xmin : xmax + 1].copy()
            ext_dq = input_model.dq[ymin : ymax + 1, xmin : xmax + 1].copy()


            #new_model = datamodels.ImageModel(data=ext_data, err=ext_err, dq=ext_dq)
            new_model = datamodels.SlitModel(data=ext_data, err=ext_err, dq=ext_dq)
            new_model.meta.wcs = subwcs
            # Not sure this makes sense for grism exposures since the trace
            # doesn't really have a footprint itself, it relates back to the
            # size of the object in the direct image. So what is really wanted
            # here?
            # util.update_s_region(new_model)
            new_model.meta.wcsinfo.spectral_order = order
            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            # nslit = obj.sid - 1  # catalog id starts at zero
            new_model.name = str(obj.sid)
            new_model.xstart = xmin + 1
            new_model.xsize = (xmax - xmin) + 1
            new_model.ystart = ymin + 1
            new_model.ysize = (ymax - ymin) + 1
            new_model.source_xpos = obj.xcentroid
            new_model.source_ypos = obj.ycentroid
            new_model.source_id = obj.sid
            new_model.bunit_data = input_model.meta.bunit_data
            new_model.bunit_err = input_model.meta.bunit_err
            slits.append(new_model)
    output_model.slits.extend(slits)
    del subwcs
    return output_model


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

