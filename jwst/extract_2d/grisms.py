#
#  Module for 2d extraction of grism spectra
#

import copy
import logging

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
    input_model : `~jwst.datamodels.CubeModel`
        The input TSO image as an instance of a CubeModel

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
    output_model : `~jwst.datamodels.SlitModel`


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
        extract_height = input_model.meta.subarray.ysize
    log.info("Setting extraction height to {}".format(extract_height))

    # Get the disperser parameters which have the wave limits
    with WavelengthrangeModel(reference_files['wavelengthrange']) as f:
        if (f.meta.instrument.name != 'NIRCAM'):
            raise ValueError("Wavelengthrange reference file not for NIRCAM!")
        if (f.meta.exposure.type != 'NRC_TSGRISM'):
            raise ValueError("Wavelengthrange reference file not for TSGRISM")
        wavelengthrange = f.wavelengthrange
        ref_extract_orders = f.extract_orders

    # this supports the user override and overriding the WFSS order extraction
    # when called from the pipeline only the first order is extracted
    if extract_orders is None:
        log.info("Using default order extraction from reference file")
        extract_orders = ref_extract_orders
        available_orders = [x[1] for x in extract_orders if x[0] == input_model.meta.instrument.filter].pop()
    else:
        if isinstance(extract_orders, list):
            if not all(isinstance(item, int) for item in extract_orders):
                raise TypeError("Expected extract_orders to be a list of ints")
        else:
            raise TypeError("Expected extract_orders to be a list of ints")
        available_orders = extract_orders

    if len(available_orders) > 1:
        raise NotImplementedError("Multiple order extraction for TSO not currently implemented")

    log.info("Extracting order: {}".format(available_orders))
    output_model = datamodels.SlitModel()
    output_model.update(input_model)
    subwcs = copy.deepcopy(input_model.meta.wcs)

    for order in available_orders:
        range_select = [(x[2], x[3]) for x in wavelengthrange if (x[0] == order and x[1] == input_model.meta.instrument.filter)]
        # All objects in the catalog will use the same filter for translation
        # that filter is the one that was used in front of the grism
        lmin, lmax = range_select.pop()

        # create the order bounding box
        source_xpos = input_model.meta.wcsinfo.crpix1 - 1  # remove fits 
        source_ypos = input_model.meta.wcsinfo.crpix2 - 1  # remove fits
        transform = input_model.meta.wcs.get_transform('full_detector', 'grism_detector')
        xmin, ymin, _ = transform(source_xpos,
                                  source_ypos,
                                  lmin,
                                  order)
        xmax, ymax, _ = transform(source_xpos,
                                  source_ypos,
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
        # in the case that it's 32 so that the height is 30 above and 34 below (in full frame)
        bump = source_ypos - 34
        extract_y_center = source_ypos - bump

        splitheight = int(extract_height / 2)
        below = extract_y_center - splitheight
        if below == 34:
            extract_y_min = 0
            extract_y_max = extract_y_center + splitheight
        elif below < 0:
            extract_y_min = 0
            extract_y_max = extract_height - 1
        else:
            extract_y_min = extract_y_center - 34  # always return source at pixel 34 in cutout.
            extract_y_max = extract_y_center + extract_height - 34 - 1

        if (extract_y_min > extract_y_max):
            raise ValueError("Something bad happened calculating extraction y-size")

        # limit the bounding boxto the detector
        ymin, ymax = (max(extract_y_min, 0), min(extract_y_max, input_model.meta.subarray.ysize))
        xmin, xmax = (max(xmin, 0), min(xmax, input_model.meta.subarray.xsize))
        log.info("xmin, xmax: {} {}  ymin, ymax: {} {}".format(xmin, xmax, ymin, ymax))

        # the order and source position are put directly into
        # the new wcs for the subarray for the forward transform
        order_model = Const1D(order)
        order_model.inverse = Const1D(order)
        tr = input_model.meta.wcs.get_transform('grism_detector', 'full_detector')
        tr = Mapping((0, 1, 0)) | Shift(xmin) & Shift(ymin) & order_model | tr
        subwcs.set_transform('grism_detector', 'full_detector', tr)

        xmin = int(xmin)
        xmax = int(xmax)
        ymin = int(ymin)
        ymax = int(ymax)

        # cut it out, keep the entire row though
        ext_data = input_model.data[:, ymin: ymax + 1, :].copy()
        ext_err = input_model.err[:, ymin: ymax + 1, :].copy()
        ext_dq = input_model.dq[:, ymin: ymax + 1, :].copy()

        log.info("WCS made explicit for order: {}".format(order))
        log.info("Trace extents are: (xmin:{}, ymin:{}), (xmax:{}, ymax:{})".format(xmin, ymin, xmax, ymax))

        if output_model.meta.model_type == "SlitModel":
            output_model.data = ext_data
            output_model.err = ext_err
            output_model.dq = ext_dq
            output_model.meta.wcs = subwcs
            output_model.meta.wcs.bounding_box = ((xmin, xmax), (ymin, ymax))
            output_model.meta.wcs.crpix2 = 34  # update for the move, vals are the same
            output_model.meta.wcsinfo.spectral_order = order
            output_model.name = str('TSO object')
            output_model.xstart = 1  # fits pixels
            output_model.xsize = ext_data.shape[2]
            output_model.ystart = 1  # fits pixels
            output_model.ysize = ext_data.shape[1]
            output_model.source_xpos = source_xpos
            output_model.source_ypos = 34
            output_model.source_id = 1
            output_model.bunit_data = input_model.meta.bunit_data
            output_model.bunit_err = input_model.meta.bunit_err

    del subwcs
    log.info("Finished extractions")
    return output_model


def extract_grism_objects(input_model,
                          grism_objects=None,
                          reference_files=None,
                          extract_orders=None,
                          use_fits_wcs=False,
                          mmag_extract=99.):
    """
    Extract 2d boxes around each objects spectra for each order.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel`
        An instance of an ImageModel, this is the grism image

    grism_objects : list(GrismObject)
        A list of GrismObjects

    reference_files : dict
        Needs to include the name of the wavelengthrange reference file


    Returns
    -------
    output_model : `~jwst.datamodels.MultiSlitModel`


    Notes
    -----
    This method supports NRC_WFSS and NIS_WFSS only

    GrismObject is a named tuple which contains distilled
    information about each catalog object. It can be created
    by calling jwst.assign_wcs.util.create_grism_bbox() which
    will return a list of GrismObjects that countain the bounding
    boxes that will be used to define the 2d extraction area.

    For each spectral order, the configuration file contains a
    magnitude-cutoff value. Sources with magnitudes fainter than the
    extraction cutoff (MMAG_EXTRACT)  will not be extracted, but are
    accounted for when computing the spectral contamination and background
    estimates. The default extraction value is 99 right now.

    The sensitivity information from the original aXe style configuration
    file needs to be modified by the passband of the filter used for
    the direct image to get the min and max wavelengths
    which correspond to t=0 and t=1, this currently has been done by the team
    and the min and max wavelengths to use to calculate t are stored in the
    grism reference file as wavelengthrange.

    Step 1: Convert the source catalog from the reference frame of the
            uberimage to that of the dispersed image.  For the Vanilla
            Pipeline we assume that the pointing information in the file
            headers is sufficient.  This will be strictly true if all images
            were obtained in a single visit (same guide stars).

    Step 2: Record source information for each object in the catalog: position
            (RA and Dec), shape (A_IMAGE, B_IMAGE, THETA_IMAGE), and all
            available magnitudes, and minimum bounding boxes

    Step 3: Compute the trace and wavelength solutions for each object in the
            catalog and for each spectral order.  Record this information.

    Step 4: Compute the WIDTH of each spectral subwindow, which may be fixed or
            variable. The cross-dispersion size is taken from the minum bounding
            box

    """
    if reference_files is None:
        raise TypeError("Expected a dictionary for reference_files")

    if grism_objects is None:
        # get the wavelengthrange reference file from the input_model
        if reference_files is None:
            raise ValueError("Need at least the dictionary of reference files")

        if (not reference_files['wavelengthrange'] or reference_files['wavelengthrange'] == 'N/A'):
            raise ValueError("Expected name of wavelengthrange reference file")
        else:
            # TODO: remove/reset the use_fits_wcs when source_catalog uses GWCS to make the catalog
            grism_objects = util.create_grism_bbox(input_model, reference_files,
                                                   extract_orders=extract_orders,
                                                   use_fits_wcs=use_fits_wcs,
                                                   mmag_extract=mmag_extract)
            log.info("Grism object list created from source catalog: {0:s}"
                     .format(input_model.meta.source_catalog.filename))

    if not isinstance(grism_objects, list):
            raise TypeError("Expected input grism objects to be a list")
    if len(grism_objects) == 0:
        raise ValueError("No grism objects created from source catalog")

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
            y, x = obj.order_bounding[order]

            # limit the boxes to the detector
            ymin = clamp(y[0], 0, input_model.meta.subarray.ysize)
            ymax = clamp(y[1], 0, input_model.meta.subarray.ysize)
            xmin = clamp(x[0], 0, input_model.meta.subarray.xsize)
            xmax = clamp(x[1], 0, input_model.meta.subarray.xsize)

            # don't extract anything that ended up with zero dimensions in one axis
            # this means that it was identified as a partial order but only on one 
            # row or column of the detector
            if (((ymax - ymin) > 0) and ((xmax - xmin) > 0)):
                # only the first two numbers in the Mapping are used
                # the order and source position are put directly into
                # the new wcs for the subarray for the forward transform
                xcenter_model = Const1D(obj.xcentroid)
                xcenter_model.inverse = Const1D(obj.xcentroid)
                ycenter_model = Const1D(obj.ycentroid)
                ycenter_model.inverse = Const1D(obj.ycentroid)
                order_model = Const1D(order)
                order_model.inverse = Const1D(order)
                tr = inwcs.get_transform('grism_detector', 'detector')
                tr = Identity(2) | Mapping((0, 1, 0, 1, 0)) | (Shift(xmin) & Shift(ymin) &
                                                               xcenter_model &
                                                               ycenter_model &
                                                               order_model) | tr

                subwcs.set_transform('grism_detector', 'detector', tr)

                log.info("Subarray extracted for obj: {} order: {}:".format(obj.sid, order))
                log.info("Subarray extents are: (xmin:{}, xmax:{}), (ymin:{}, ymax:{})".format(xmin, xmax, ymin, ymax))

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
                log.debug("sid {} xcen {}  ycen {}".format(obj.sid, obj.xcentroid, obj.ycentroid))
                new_model.name = str(obj.sid)
                new_model.xstart = xmin + 1
                new_model.xsize = (xmax - xmin)
                new_model.ystart = ymin + 1
                new_model.ysize = (ymax - ymin)
                new_model.source_xpos = float(obj.xcentroid)
                new_model.source_ypos = float(obj.ycentroid)
                new_model.source_id = obj.sid
                new_model.bunit_data = input_model.meta.bunit_data
                new_model.bunit_err = input_model.meta.bunit_err
                new_model.meta.wcs.bounding_box = ((xmin, xmax), (ymin, ymax))
                slits.append(new_model)

    output_model.slits.extend(slits)
    del subwcs
    return output_model


def clamp(value, minval, maxval):
    """
    Return the value clipped between minval and maxval.

    Parameters
    ----------
    value: float
        The value to limit
    minval: float
        The minimal acceptable value
    maxval: float
        The maximum acceptable value

    Returns
    -------
    value: float
        The value that falls within the min-max range or the minimum limit
    """
    return max(minval, min(value, maxval))


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
