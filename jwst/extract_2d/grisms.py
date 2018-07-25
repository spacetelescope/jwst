#
#  Module for 2d extraction of grism spectra
#

import logging
import copy
from astropy.modeling.models import Shift, Const1D, Mapping, Identity

from .. import datamodels
from ..assign_wcs import util
from ..datamodels import CubeModel, WavelengthrangeModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract_tso_object(input_model,
                       reference_files=None,
                       extract_height=None,
                       extract_orders=None):
    """
    Extract the spectrum for a NIRCAM TSO observation.

    Parameters
    ----------
    input_model : jwst.datamodels.CubeModel

    reference_files : dict
        Needs to include the name of the wavelengthrange reference file

    extract_height: int
        The extraction height, in total, for the spectra in the
        cross-dispersion direction. If this is other than None,
        it will override the team default of 64 pixels. The team
        wants the source centered near row 34, so the extraction
        height is not the same on either size of the central pixel.

    extract_orders: list[ints]
        This is an optional parameter that will override the
        orders specified for extraction in the wavelengthrange
        reference file.

    Returns
    -------
    output_model : jwst.datamodels.MultiSlitModel()


    Notes
    -----
    This method supports NRC_TSGRISM only, where only one
    bright object is considered in the field, so there's
    no catalog of sources and the object is assumed to
    have been observed at the aperture location crpix1/2.

    GrismObject is a named tuple which contains distilled
    information about each catalog object and the bounding
    boxes that will be used to define the 2d extraction area.

    Since this mode has a single known source location the utilities
    used in the WFSS modes are overkill, instead, similar structures
    are created during the extrac 2d process and then directly used.

    https://jwst-docs.stsci.edu/display/JTI/NIRCam+Grism+Time+Series
    """

    if not isinstance(reference_files, dict):
        raise TypeError("Expected a dictionary for reference_files")
    if 'wavelengthrange' not in reference_files.keys():
        raise KeyError("No wavelengthrange reference file specified")

    if not isinstance(input_model, CubeModel):
        raise TypeError('The input data model is not a CubeModel.')

    if extract_height is None:
        extract_height = 64  # set by the teams

    log.info("Extracting into a MultiSlitModel")
    output_model = datamodels.MultiSlitModel()
    output_model.update(input_model)

    subwcs = copy.deepcopy(input_model.meta.wcs)

    # Get the disperser parameters which have the wave limits
    with WavelengthrangeModel(reference_files['wavelengthrange']) as f:
        if (f.meta.instrument.name != 'NIRCAM'):
            raise ValueError("Wavelengthrange reference file not for NIRCAM!")
        wavelengthrange = f.wavelengthrange
        ref_extract_orders = f.extract_orders

    # this supports the user override
    if extract_orders is None:
        extract_orders = ref_extract_orders

    available_orders = [x[1] for x in extract_orders if x[0] == input_model.meta.instrument.filter].pop()
    length, _, _ = input_model.data.shape

    slits = []
    for trace_image in range(length):
        for order in available_orders:
            range_select = [(x[2], x[3]) for x in wavelengthrange if (x[0] == order and x[1] == input_model.meta.instrument.filter)]
            # All objects in the catalog will use the same filter for translation
            # that filter is the one that was used in front of the grism
            lmin, lmax = range_select.pop()

            # create the order bounding box
            transform = input_model.meta.wcs.get_transform('full_detector', 'grism_detector')
            xmin, ymin, _ = transform(input_model.meta.wcsinfo.crpix1,
                                      input_model.meta.wcsinfo.crpix2,
                                      lmin,
                                      order)
            xmax, ymax, _ = transform(input_model.meta.wcsinfo.crpix1,
                                      input_model.meta.wcsinfo.crpix2,
                                      lmax,
                                      order)

            # Add the shift to the lower corner to each subarray WCS object
            # The shift should just be the lower bounding box corner
            # also replace the object center location inputs to the GrismDispersion
            # model with the known object center and order information (in pixels of direct image)
            # This is changes the user input to the model from (x,y,x0,y0,order) -> (x,y)
            #
            # The bounding boxes here are also limited to the size of the detector in the
            # dispersion direction and 64 pixels in the cross-dispersion 
            # The check for boxes entirely off the detector is done in create_grism_bbox right now

            # The team wants the object to fall near  row 34 for all cutouts,
            # but the default cutout height is 64pixel (32 on either side)
            # so use crpix2 when it equals 34, but  bump the ycenter by 2 pixel
            # in other cases  so that the height is 30 above and 34 below (in full frame)
            extract_y_center = input_model.meta.wcsinfo.crpix2
            if input_model.meta.wcsinfo.crpix2 != 34:
                extract_y_center = extract_y_center - 2
            extract_y_min = extract_y_center - 34
            extract_y_max = extract_y_center + 30

            # limit the boxes to the detector
            ymin, ymax = (max(extract_y_min, 0), min(extract_y_max, input_model.meta.subarray.ysize))
            xmin, xmax = (max(xmin, 0), min(xmax, input_model.meta.subarray.xsize))
            log.info("xmin, xmax: {} {}  ymin, ymax: {} {}".format(xmin, xmax, ymin, ymax))
            # only the first two numbers in the Mapping are used
            # the order and source position are put directly into
            # the new wcs for the subarray for the forward transform
            order_model = Const1D(order)
            order_model.inverse = Const1D(order)
            tr = input_model.meta.wcs.get_transform('grism_detector', 'full_detector')
            tr = Mapping((0, 1, 0)) | (Shift(xmin) & Shift(ymin) & order_model) | tr
            subwcs.set_transform('grism_detector', 'full_detector', tr)
            log.info("Extracting from image: {}".format(trace_image))
            log.info("Subarray extracted for order: {}:".format(order))
            log.info("Subarray extents are: (xmin:{}, ymin:{}), (xmax:{}, ymax:{})".format(xmin, ymin, xmax, ymax))

            xmin = int(xmin)
            xmax = int(xmax)
            ymin = int(ymin)
            ymax = int(ymax)

            # cut it out
            ext_data = input_model.data[trace_image][ymin: ymax + 1, xmin: xmax + 1].copy()
            ext_err = input_model.err[trace_image][ymin: ymax + 1, xmin: xmax + 1].copy()
            ext_dq = input_model.dq[trace_image][ymin: ymax + 1, xmin: xmax + 1].copy()
            new_model = datamodels.SlitModel(data=ext_data, err=ext_err, dq=ext_dq)
            new_model.meta.wcs = subwcs

            new_model.meta.wcsinfo.spectral_order = order
            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            # nslit = obj.sid - 1  # catalog id starts at zero
            new_model.name = str(1)
            new_model.xstart = xmin + 1
            new_model.xsize = (xmax - xmin) + 1
            new_model.ystart = ymin + 1
            new_model.ysize = (ymax - ymin) + 1
            new_model.source_xpos = input_model.meta.wcsinfo.crpix1
            new_model.source_ypos = input_model.meta.wcsinfo.crpix2
            new_model.source_id = 1
            new_model.bunit_data = input_model.meta.bunit_data
            new_model.bunit_err = input_model.meta.bunit_err
            slits.append(new_model)

    output_model.slits.extend(slits)
    del subwcs
    log.info("Finished extractions")
    return output_model


def extract_grism_objects(input_model,
                          grism_objects=None,
                          reference_files=None):
    """
    Extract 2d boxes around each objects spectra for each order.

    Parameters
    ----------
    input_model : jwst.datamodels.ImageModel

    grism_objects : list(GrismObject)
        A list of GrismObjects

    reference_files : dict
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
    if reference_files is None:
        raise TypeError("Expected a dictionary for reference_files")

    if grism_objects is None:
        # get the wavelengthrange reference file from the input_model
        if (not reference_files['wavelengthrange'] or reference_files['wavelengthrange'] == 'N/A'):
            raise ValueError("Expected name of wavelengthrange reference file")
        else:
            grism_objects = util.create_grism_bbox(input_model, reference_files)
            log.info("Grism object list created from source catalog: {0:s}"
                     .format(input_model.meta.source_catalog.filename))

    if not isinstance(grism_objects, list):
            raise ValueError("Expected input grism objects to be a list")

    log.info("Extracting grism objects into MultiSlitModel")
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
            log.info("Subarray extents are: (xmin:{}, ymin:{}), (xmax:{}, ymax:{})".format(xmin, ymin, xmax, ymax))

            ext_data = input_model.data[ymin: ymax + 1, xmin: xmax + 1].copy()
            ext_err = input_model.err[ymin: ymax + 1, xmin: xmax + 1].copy()
            ext_dq = input_model.dq[ymin: ymax + 1, xmin: xmax + 1].copy()

            # new_model = datamodels.ImageModel(data=ext_data, err=ext_err, dq=ext_dq)
            new_model = datamodels.SlitModel(data=ext_data, err=ext_err, dq=ext_dq)
            new_model.meta.wcs = subwcs
            new_model.meta.wcsinfo.spectral_order = order

            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            # nslit = obj.sid - 1  # catalog id starts at zero
            new_model.name = str(obj.sid)
            new_model.xstart = xmin + 1
            new_model.xsize = (xmax - xmin)
            new_model.ystart = ymin + 1
            new_model.ysize = (ymax - ymin)
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

