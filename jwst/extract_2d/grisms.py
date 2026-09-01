"""Functions for 2D extraction of grism spectra."""

import copy
import logging

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.modeling import bind_bounding_box
from astropy.modeling.models import Const1D, Mapping, Shift
from gwcs.utils import to_index
from gwcs.wcstools import grid_from_bounding_box
from stcal.alignment.util import wcs_bbox_from_shape
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import WavelengthrangeModel

from jwst.assign_wcs import util
from jwst.lib.catalog_utils import read_source_catalog
from jwst.lib.stripe_utils import generate_substripe_ranges

log = logging.getLogger(__name__)

__all__ = [
    "extract_tso_object",
    "extract_grism_objects",
    "compute_dispersion",
    "compute_tso_wavelength_array",
    "compute_wfss_wavelength",
]


def build_grism_submodel(
    sub_model,
    input_model,
    xmin,
    xmax,
    ymin,
    ymax,
    subwcs,
    compute_wavelength,
    order,
    name="1",
    source_xpos=None,
    source_ypos=None,
):
    """
    Build a grism model from the input data.

    Parameters
    ----------
    sub_model : `~stdatamodels.jwst.datamodels.SlitModel`
        The data model to be filled with arrays and WCS information.
    input_model : `~stdatamodels.jwst.datamodels.CubeModel`
        The parent model from which the 2D extraction is taken.
    xmin : int
        The minimum x pixel column value for the extracted region.
    xmax : int
        The maximum x pixel column value for the extracted region.
    ymin : int
        The minimum y pixel column value for the extracted region.
    ymax : int
        The maximum y pixel column value for the extracted region.
    subwcs : `~gwcs.wcs.WCS`
        The WCS object from the parent model, modified to fit the
        extracted region.
    compute_wavelength : bool
        If `True`, compute the wavelength array of the extracted region.
    order : int
        The spectral order of the extracted region.
    name : str, optional
        The name of the extracted region; typically a placeholder
        for NRC_TSGRISM data but will be the stripe number for DHS.
    source_xpos : float, optional
        The x position of the source in the direct image frame (0-indexed).
        When provided, sets ``sub_model.source_xpos`` and updates
        ``sub_model.meta.wcsinfo.siaf_xref_sci``.
    source_ypos : float, optional
        The y position of the source in the direct image frame (0-indexed).
        When provided, sets ``sub_model.source_ypos``.
    """
    # Cut out the subarray from the input data arrays
    ext_data = input_model.data[..., ymin : ymax + 1, xmin : xmax + 1].copy()
    ext_err = input_model.err[..., ymin : ymax + 1, xmin : xmax + 1].copy()
    ext_dq = input_model.dq[..., ymin : ymax + 1, xmin : xmax + 1].copy()
    if input_model.var_poisson is not None and np.size(input_model.var_poisson) > 0:
        var_poisson = input_model.var_poisson[..., ymin : ymax + 1, xmin : xmax + 1].copy()
    else:
        var_poisson = None
    if input_model.var_rnoise is not None and np.size(input_model.var_rnoise) > 0:
        var_rnoise = input_model.var_rnoise[..., ymin : ymax + 1, xmin : xmax + 1].copy()
    else:
        var_rnoise = None
    if input_model.var_flat is not None and np.size(input_model.var_flat) > 0:
        var_flat = input_model.var_flat[..., ymin : ymax + 1, xmin : xmax + 1].copy()
    else:
        var_flat = None

    # Finish populating the output model and meta data
    sub_model.data = ext_data
    sub_model.err = ext_err
    sub_model.dq = ext_dq
    sub_model.var_poisson = var_poisson
    sub_model.var_rnoise = var_rnoise
    sub_model.var_flat = var_flat
    sub_model.meta.wcs = subwcs
    sub_model.meta.wcs.bounding_box = wcs_bbox_from_shape(ext_data.shape)
    if compute_wavelength:
        sub_model.wavelength = compute_tso_wavelength_array(sub_model)
    sub_model.meta.wcsinfo.spectral_order = order
    sub_model.meta.wcsinfo.dispersion_direction = input_model.meta.wcsinfo.dispersion_direction
    sub_model.meta.instrument.name = "NIRCAM"
    sub_model.name = name
    sub_model.source_type = input_model.meta.target.source_type
    sub_model.source_name = input_model.meta.target.catalog_name
    sub_model.source_alias = input_model.meta.target.proposer_name
    sub_model.xstart = 1  # FITS pixels are 1-indexed
    sub_model.xsize = ext_data.shape[-1]
    sub_model.ystart = ymin + 1  # FITS pixels are 1-indexed
    sub_model.ysize = ext_data.shape[-2]
    if source_xpos is not None:
        sub_model.source_xpos = source_xpos
    if source_ypos is not None:
        sub_model.source_ypos = source_ypos

        ra, dec, _, _ = subwcs(ext_data.shape[-1] / 2, source_ypos)
        sub_model.source_ra = ra
        sub_model.source_dec = dec

    sub_model.source_id = 1
    sub_model.meta.bunit_data = input_model.meta.bunit_data
    sub_model.meta.bunit_err = input_model.meta.bunit_err
    if getattr(input_model, "int_times", None) is not None:
        sub_model.int_times = input_model.int_times.copy()


def _set_tso_subwcs_transform(input_model, subwcs, xstart, ymin, order):
    """Make grism to direct image transform for the subwcs."""  # numpydoc ignore:RT01
    order_model = Const1D(order)
    order_model.inverse = Const1D(order)
    tr = input_model.meta.wcs.get_transform("grism_detector", "direct_image")
    tr = Mapping((0, 1, 0)) | Shift(xstart) & Shift(ymin) & order_model | tr
    subwcs.set_transform("grism_detector", "direct_image", tr)


def extract_tso_object(
    input_model,
    reference_files=None,
    tsgrism_extract_height=None,
    extract_orders=None,
    compute_wavelength=True,
):
    """
    Extract the spectrum for a NIRCam TSO grism observation.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.CubeModel` or \
                  `~stdatamodels.jwst.datamodels.ImageModel`
        The input TSO data can be a cube (3D) or an image (2D).

    reference_files : dict
        This dictionary must contain the name of the
        WAVELENGTHRANGE reference file.

    tsgrism_extract_height : int, optional
        The extraction height, in total, for the spectrum in the
        cross-dispersion direction. If this is other than None,
        it will override the default of 64 pixels. The instrument team
        wants the source centered near row 34, so the extraction
        height is not the same on either size of the central row.

    extract_orders : list of int, optional
        Overrides the orders specified for extraction in the
        WAVELENGTHRANGE reference file.

    compute_wavelength : bool, optional
        Compute a wavelength array for the datamodel.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.SlitModel`
        Output model containing extracted spectrum.

    Notes
    -----
    This method supports NRC_TSGRISM only, where only one bright object is
    considered in the field, so there's no catalog of sources and the object
    is assumed to have been observed at the aperture reference position.
    The aperture reference location is read during level-1b (uncal) product
    creation by the "set_telescope_pointing" script from the SIAF entries
    XSciRef and YSciRef (reference location in the science frame) and saved as
    "meta.wcsinfo.siaf_xref_sci" and "meta.wcsinfo.siaf_yref_sci" (FITS header
    keywords XREF_SCI and YREF_SCI).

    Because this mode has a single known source location, the utilities used
    in the WFSS modes are overkill. Instead, similar structures are created
    during the extract2d process and then directly used.

    For more information on the NRC_TSGRISM mode, see:
    https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-observing-modes/nircam-time-series-observations/nircam-grism-time-series
    """
    # Check for reference files
    if not isinstance(reference_files, dict):
        raise TypeError("Expected a dictionary for reference_files")

    # Check for wavelengthrange reference file
    if "wavelengthrange" not in reference_files:
        raise KeyError("No wavelengthrange reference file specified")

    # Get the disperser parameters that have the wave limits
    with WavelengthrangeModel(reference_files["wavelengthrange"]) as f:
        if f.meta.instrument.name != "NIRCAM" or f.meta.exposure.type != "NRC_TSGRISM":
            raise ValueError("Wavelengthrange reference file is not for NIRCAM TSGRISM mode!")
        wavelengthrange = f.wavelengthrange
        ref_extract_orders = f.extract_orders

    # If user-supplied spectral orders are not provided,
    # default to extracting only the 1st order
    if extract_orders is None:
        log.info("Using default order extraction from reference file")
        extract_orders = ref_extract_orders
        available_orders = [
            x[1] for x in extract_orders if x[0] == input_model.meta.instrument.filter
        ].pop()
    else:
        if not isinstance(extract_orders, list) or not all(
            isinstance(item, int) for item in extract_orders
        ):
            raise TypeError("Expected extract_orders to be a list of integers.")
        available_orders = extract_orders

    if len(available_orders) > 1:
        raise NotImplementedError("Multiple order extraction for TSO is not currently implemented.")

    # Check for the existence of the aperture reference location meta data
    if (
        input_model.meta.wcsinfo.siaf_xref_sci is None
        or input_model.meta.wcsinfo.siaf_yref_sci is None
    ):
        raise ValueError("XREF_SCI and YREF_SCI are required for TSO mode.")

    # Split the logic on DHS vs. non-DHS data
    if "DHS" in input_model.meta.subarray.name.upper():
        output_model = _extract_tso_dhs_object(
            input_model, wavelengthrange, available_orders, compute_wavelength
        )
    else:
        output_model = _extract_tso_tsgrism_object(
            input_model,
            wavelengthrange,
            available_orders,
            compute_wavelength,
            tsgrism_extract_height=tsgrism_extract_height,
        )
    log.info("Finished extraction")
    return output_model


def _extract_tso_tsgrism_object(
    input_model, wavelengthrange, available_orders, compute_wavelength, tsgrism_extract_height=None
):
    # Processing non-DHS NRC_TSGRISM data
    # Create the extracted output as a SlitModel
    log.info(f"Extracting order: {available_orders}")
    output_model = datamodels.SlitModel()
    output_model.update(input_model)

    subwcs = copy.deepcopy(input_model.meta.wcs)

    # TODO: Moved this from above to segment DHS vs. non-DHS data - check compat
    # If an extraction height is not supplied, default to entire
    # cross-dispersion size of the data array
    if tsgrism_extract_height is None:
        tsgrism_extract_height = input_model.meta.subarray.ysize
    log.info(f"Setting extraction height to {tsgrism_extract_height}")

    # Loop over spectral orders
    for order in available_orders:
        _, _, _, lmin, lmax = [
            x
            for x in wavelengthrange
            if (x[0] == order and x[2] == input_model.meta.instrument.filter)
        ][0]

        # Source xpos can be taken as the reference x position
        source_xpos = input_model.meta.wcsinfo.siaf_xref_sci - 1

        # Take the input source position to be row 34, always
        input_ypos = 34.0

        if tsgrism_extract_height >= 64:
            # For the default 64-pixel cutout or larger, start at 0 so
            # that source position stays at row 34
            extract_y_min = 0
        else:
            # Otherwise, center the extraction on the source position
            # This may require a custom extract1d reference file to center
            # the spectrum correctly.
            extract_y_min = input_ypos - tsgrism_extract_height / 2
        extract_y_max = extract_y_min + tsgrism_extract_height - 1

        # Limit the bounding box to the detector edges
        # Note: in contrast to WFSS modes, the instrument team requested that the
        # entire detector be extracted in the x-direction, rather than determining
        # the min/max values for x from the min/max wavelength values in the WCS.
        ny, nx = input_model.data.shape[-2:]
        ymin = int(np.clip(extract_y_min, 0, ny - 1))
        ymax = int(np.clip(extract_y_max, ymin, ny - 1))
        xmin = 0
        xmax = nx - 1

        # Output source position may be shifted
        source_ypos = input_ypos - ymin

        # The order and source position are put directly into the new WCS of the subarray
        # for the forward transform.
        _set_tso_subwcs_transform(input_model, subwcs, xmin, ymin, order)

        log.info(f"WCS made explicit for order: {order}")
        log.info(f"Extraction limits: (xmin: {xmin}, ymin: {ymin}), (xmax: {xmax}, ymax: {ymax})")

        build_grism_submodel(
            output_model,
            input_model,
            xmin,
            xmax,
            ymin,
            ymax,
            subwcs,
            compute_wavelength,
            order,
            name="1",
            source_xpos=source_xpos,
            source_ypos=source_ypos,
        )
    del subwcs
    return output_model


def _extract_tso_dhs_object(
    input_model,
    wavelengthrange,
    available_orders,
    compute_wavelength=True,
):
    """
    Extract the spectra for a NIRCam DHS TSO observation.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.CubeModel` or \
                  `~stdatamodels.jwst.datamodels.ImageModel`
        The input TSO DHS data.
    wavelengthrange : list
        The wavelength range table from the WAVELENGTHRANGE reference file.
    available_orders : list of int
        The spectral orders to extract.
    compute_wavelength : bool
        Compute a wavelength array for the datamodel.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        Output model with one slit per stripe.
    """
    output_model = datamodels.MultiSlitModel()
    output_model.update(input_model)

    data_shape = input_model.data.shape
    xx, yy = np.meshgrid(np.arange(data_shape[-1]), np.arange(data_shape[-2]))
    fwd_xfrm = input_model.meta.wcs.get_transform("grism_detector", "direct_image")
    all_stripes = fwd_xfrm(xx, yy, np.ones_like(xx))[-1]
    if "LONG" in input_model.meta.instrument.detector.upper():
        # Because nrcalong DHS repeats reads of the same detector position
        # for all stripes, generate a list of stripe numbers from the subarray
        # name rather than unique regions values.
        subarray_stripenum = int(input_model.meta.subarray.name.split("STRIPE")[1][0])
        stripe_set = np.array(range(subarray_stripenum)) + 1
        sub_ranges = generate_substripe_ranges(input_model, science_frame=True)["subarray"]
    else:
        # For short wavelength detectors, use region values directly
        stripe_set = np.unique(all_stripes[~np.isnan(all_stripes)].astype(int))

    for i, stripe_id in enumerate(stripe_set):
        for order in available_orders:
            sub_model = datamodels.SlitModel()
            subwcs = copy.deepcopy(input_model.meta.wcs)

            fieldpoint_idx = 1
            filter_idx = 2

            waverange_match = [
                x
                for x in wavelengthrange
                if (x[0] == order and x[filter_idx] == input_model.meta.instrument.filter)
            ]

            if len(waverange_match) > 1:
                _, _, _, lmin, lmax = [
                    x
                    for x in waverange_match
                    if x[fieldpoint_idx] in input_model.meta.aperture.pps_name
                ][0]
            else:
                _, _, _, lmin, lmax = waverange_match[0]

            # Find extent of stripe as defined by regions
            # For nrcalong, the regions is not helpful, so rely on
            # ranges generated by readout recreation
            if "LONG" in input_model.meta.instrument.detector.upper():
                stripe_x = xx
                # Range generated from array slice, which causes unwanted
                # extra row - drop it here.
                stripe_y = np.array([sub_ranges[i][0], sub_ranges[i][1] - 1])
            else:
                stripe_x = np.where(all_stripes == stripe_id, xx, np.nan)
                stripe_y = np.where(all_stripes == stripe_id, yy, np.nan)
            stripe_xmin = np.nanmin(stripe_x)
            stripe_xmax = np.nanmax(stripe_x)
            stripe_ymin = np.nanmin(stripe_y)
            stripe_ymax = np.nanmax(stripe_y)

            xmin, xmax = (
                max(stripe_xmin, 0),
                min(stripe_xmax, input_model.meta.subarray.xsize),
            )
            ymin, ymax = (
                max(stripe_ymin, 0),
                min(stripe_ymax, input_model.meta.subarray.ysize),
            )

            _set_tso_subwcs_transform(input_model, subwcs, xmin, ymin, order)

            xmin = int(xmin)
            xmax = int(xmax)
            ymin = int(ymin)
            ymax = int(ymax)

            log.info(f"WCS made explicit for stripe {stripe_id}, order {order}.")
            log.info(
                f"Extraction limits: (xmin: {xmin}, ymin: {ymin}), (xmax: {xmax}, ymax: {ymax})"
            )

            build_grism_submodel(
                sub_model,
                input_model,
                xmin,
                xmax,
                ymin,
                ymax,
                subwcs,
                compute_wavelength,
                order,
                name=str(stripe_id),
            )
            output_model.slits.append(sub_model)
    if getattr(input_model, "int_times", None) is not None:
        output_model.int_times = input_model.int_times.copy()
    return output_model


def extract_grism_objects(
    input_model,
    grism_objects=None,
    reference_files=None,
    extract_orders=None,
    source_ids=None,
    source_ra=None,
    source_dec=None,
    max_sep=None,
    mmag_extract=None,
    compute_wavelength=True,
    wfss_extract_half_height=None,
    nbright=None,
):
    """
    Extract 2D boxes around each objects spectra for each order.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.ImageModel`
        Model of the grism image.

    grism_objects : list of `~stdatamodels.jwst.transforms.GrismObject`
        A list of grism objects.

    reference_files : dict
        This dictionary must contain the name of the
        WAVELENGTHRANGE reference file.

    extract_orders : int
        Spectral orders to extract.

    source_ids : list
        List of source IDs to extract.

    source_ra : list of float
        Source right ascensions to be processed. The nearest matching source to each RA/Dec pair
        will be extracted. If both ``source_ids`` and ``source_ra``/``source_dec`` are provided,
        the lists will be combined and their union extracted.

    source_dec : list of float
        Source declinations to be processed, must have same length as ``source_ra``.

    max_sep : float
        Radius in arcseconds within which ``source_ra`` and ``source_dec`` will be matched
        to sources in the catalog. If no source is found within this radius, a warning
        will be emitted and no source will be extracted corresponding to that RA, Dec pair.

    mmag_extract : float
        The minimum magnitude extraction cutoff. Sources fainter than this
        will not be extracted.

    compute_wavelength : bool
        Compute a wavelength array for the datamodel.

    wfss_extract_half_height : int
        Cross-dispersion extraction half height in pixels.
        Overwrites the computed extraction height.

    nbright : int
        Number of brightest objects to extract.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        Output model of extracted spectra.

    Notes
    -----
    This method supports WFSS modes only.

    `~stdatamodels.jwst.transforms.GrismObject`
    is a named tuple which contains distilled
    information about each catalog object. It can be created
    by calling :func:`~jwst.assign_wcs.util.create_grism_bbox` which
    will return a list of `~stdatamodels.jwst.transforms.GrismObject`
    that contains the bounding
    boxes that will be used to define the 2D extraction area.

    For each spectral order, the configuration file contains a
    magnitude-cutoff value. The total list of objects to extract is limited
    by both MMAG_EXTRACT and NBRIGHT. Sources with magnitudes fainter than the
    extraction cutoff (MMAG_EXTRACT) will not be extracted, but are
    accounted for when computing the spectral contamination and background
    estimates; the default value is 99.
    NBRIGHT further limits the list to the NBRIGHT brightest objects;
    the default value is 999.

    The sensitivity information from the original aXe-style configuration
    file needs to be modified by the passband of the filter used for
    the direct image to get the min and max wavelengths
    which correspond to ``t=0`` and ``t=1``.
    The min and max wavelengths used to calculate ``t`` are stored in the
    grism WAVELENGTHRANGE reference file.

    1. Convert the source catalog from the reference frame of the
       uber-image to that of the dispersed image.
       We assume that the pointing information in the file
       headers is sufficient.  This will be strictly true if all images
       were obtained in a single visit (same guide stars).

    2. Record source information for each object in the catalog: position
       (RA, Dec), shape (A_IMAGE, B_IMAGE, THETA_IMAGE), and all
       available magnitudes, and minimum bounding boxes.

    3. Compute the trace and wavelength solutions for each object in the
       catalog and for each spectral order.  Record this information.

    4. Compute the WIDTH of each spectral subwindow, which may be fixed or
       variable. The cross-dispersion size is taken from the minimum
       bounding box.

    Each of the virtual slits in the output
    `~stdatamodels.jwst.datamodels.MultiSlitModel` will have its own
    WCS object that is a copy of the input model's WCS, but with an additional
    transform from "grism_slit" to "grism_detector" prepended to it; this
    transform encodes a shift to the center of the slit and a binding to the
    slit's bounding box.
    """
    if reference_files is None or not reference_files:
        raise TypeError("Expected a dictionary for reference_files")

    if grism_objects is None:
        # get the wavelengthrange reference file from the input_model
        if "wavelengthrange" not in reference_files or reference_files["wavelengthrange"] in [
            "N/A",
            "",
        ]:
            raise ValueError("Expected name of wavelengthrange reference file")

        source_ids = radec_to_source_ids(
            input_model.meta.source_catalog, source_ids, source_ra, source_dec, max_sep=max_sep
        )
        grism_objects = util.create_grism_bbox(
            input_model,
            reference_files,
            extract_orders=extract_orders,
            source_ids=source_ids,
            mmag_extract=mmag_extract,
            wfss_extract_half_height=wfss_extract_half_height,
            nbright=nbright,
        )
        log.info(
            f"Grism object list created from source catalog: {input_model.meta.source_catalog}"
        )

    if not isinstance(grism_objects, list):
        raise TypeError("Expected input grism objects to be a list")
    if len(grism_objects) == 0:
        raise ValueError("No grism objects created from source catalog")

    log.info(f"Extracting {len(grism_objects)} grism objects")
    output_model = datamodels.MultiSlitModel(validate_on_assignment=False)
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
            ymin = np.clip(y[0], 0, input_model.meta.subarray.ysize)
            log.debug(f"YYY, {y}, {ymin}")
            ymax = np.clip(y[1], 0, input_model.meta.subarray.ysize)
            xmin = np.clip(x[0], 0, input_model.meta.subarray.xsize)
            xmax = np.clip(x[1], 0, input_model.meta.subarray.xsize)

            # don't extract anything that ended up with zero dimensions in one axis
            # this means that it was identified as a partial order but only on one
            # row or column of the detector
            if ymax - ymin > 0 and xmax - xmin > 0:
                subwcs = copy.deepcopy(inwcs)
                log.info(f"Subarray extracted for obj: {obj.sid} order: {order}:")
                log.info(
                    f"Subarray extents are: (xmin:{xmin}, xmax:{xmax}), (ymin:{ymin}, ymax:{ymax})"
                )

                # only the first two numbers in the Mapping are used
                # the order and source position are put directly into
                # the new wcs for the subarray for the forward transform
                xcenter_model = Const1D(obj.xcentroid)
                xcenter_model.inverse = Const1D(obj.xcentroid)

                ycenter_model = Const1D(obj.ycentroid)
                ycenter_model.inverse = Const1D(obj.ycentroid)

                order_model = Const1D(order)
                order_model.inverse = Const1D(order)

                y_slice = slice(to_index(ymin), to_index(ymax) + 1)
                x_slice = slice(to_index(xmin), to_index(xmax) + 1)

                ext_data = input_model.data[y_slice, x_slice].copy()
                ext_err = input_model.err[y_slice, x_slice].copy()
                ext_dq = input_model.dq[y_slice, x_slice].copy()
                if input_model.var_poisson is not None and np.size(input_model.var_poisson) > 0:
                    var_poisson = input_model.var_poisson[y_slice, x_slice].copy()
                else:
                    var_poisson = None
                if input_model.var_rnoise is not None and np.size(input_model.var_rnoise) > 0:
                    var_rnoise = input_model.var_rnoise[y_slice, x_slice].copy()
                else:
                    var_rnoise = None
                if input_model.var_flat is not None and np.size(input_model.var_flat) > 0:
                    var_flat = input_model.var_flat[y_slice, x_slice].copy()
                else:
                    var_flat = None

                # Add a new transform to the WCS that shifts to the center of the virtual slit
                # This needs to be separated from the "grism_detector"/("dispersed_detector")
                # to "detector" transform  because the un-shifted "grism_detector" to "detector"
                # transform is used by wfss_contam

                tr = Mapping((0, 1, 0, 0, 0)) | (
                    Shift(xmin) & Shift(ymin) & xcenter_model & ycenter_model & order_model
                )
                bind_bounding_box(
                    tr, util.transform_bbox_from_shape(ext_data.shape, order="F"), order="F"
                )

                grism_slit = copy.deepcopy(subwcs.grism_detector)
                grism_slit.name = "grism_slit"
                subwcs.insert_frame(
                    input_frame=grism_slit, output_frame="grism_detector", transform=tr
                )
                # Force the pipelines to share their grism_detector-world transforms.
                # We want that transform to be serialized just once on save
                # instead of copied a bunch of times.
                # It was found that validation of all those copies is very slow, and inflates
                # the file size unnecessarily.
                subwcs.pipeline[1:] = inwcs.pipeline[:]

                new_slit = datamodels.SlitModel(
                    data=ext_data,
                    err=ext_err,
                    dq=ext_dq,
                    var_poisson=var_poisson,
                    var_rnoise=var_rnoise,
                    var_flat=var_flat,
                    validate_on_assignment=False,  # for runtime
                )

                new_slit.meta.wcsinfo.spectral_order = order
                new_slit.meta.wcsinfo.dispersion_direction = (
                    input_model.meta.wcsinfo.dispersion_direction
                )
                new_slit.meta.wcsinfo.specsys = input_model.meta.wcsinfo.specsys
                new_slit.meta.coordinates = input_model.meta.coordinates
                new_slit.meta.wcs = subwcs

                if compute_wavelength:
                    log.debug("Computing wavelengths")
                    new_slit.wavelength = compute_wfss_wavelength(new_slit)

                # set x/ystart values relative to the image (screen) frame.
                # The overall subarray offset is recorded in model.meta.subarray.
                # nslit = obj.sid - 1  # catalog id starts at zero
                new_slit.name = f"{obj.sid}"
                new_slit.is_extended = obj.is_extended
                new_slit.xstart = to_index(xmin) + 1  # fits pixels
                new_slit.xsize = ext_data.shape[1]
                new_slit.ystart = to_index(ymin) + 1  # fits pixels
                new_slit.ysize = ext_data.shape[0]
                new_slit.source_xpos = float(obj.xcentroid)
                new_slit.source_ypos = float(obj.ycentroid)
                new_slit.source_id = obj.sid
                new_slit.source_dec = obj.sky_centroid.dec.value
                new_slit.source_ra = obj.sky_centroid.ra.value
                new_slit.meta.bunit_data = input_model.meta.bunit_data
                new_slit.meta.bunit_err = input_model.meta.bunit_err
                slits.append(new_slit)
    output_model.slits.extend(slits)

    # update s_region of 0th slit to match input model
    if output_model.slits:
        output_model.slits[0].meta.wcsinfo.s_region = input_model.meta.wcsinfo.s_region

    # In the case that there are no spectra to extract deleting the variables
    # will fail so add the try block.
    try:
        del subwcs
    except UnboundLocalError:
        pass
    try:
        del new_slit
    except UnboundLocalError:
        pass
    # del subwcs
    # del new_slit
    log.info("Finished extractions")
    return output_model


def compute_dispersion(wcs):
    """
    Compute the pixel dispersion.

    Make a model for the pixel dispersion from the ``grismconf`` specs.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.

    Returns
    -------
    dispersion : ndarray
        The pixel dispersion in meters.
    """
    raise NotImplementedError


def compute_tso_wavelength_array(slit):
    """
    Compute the wavelength array for a slit with WCS.

    Parameters
    ----------
    slit : `~stdatamodels.jwst.datamodels.SlitModel`
        JWST slit datamodel containing a ``meta.wcs`` that is a
        `~gwcs.wcs.WCS` object

    Returns
    -------
    wavelength : ndarray
        The wavelength array
    """
    wcs = slit.meta.wcs
    full_transform = slit.meta.wcs.forward_transform
    x, y = grid_from_bounding_box(wcs.bounding_box)
    wavelength = full_transform(x, y)[2]
    return wavelength


def compute_wfss_wavelength(slit):
    """
    Compute the wavelength array for a slit with WCS.

    Parameters
    ----------
    slit : `~stdatamodels.jwst.datamodels.SlitModel`
        JWST slit datamodel containing a ``meta.wcs`` that is a
        `~gwcs.wcs.WCS` object

    Returns
    -------
    wavelength : ndarray
        The wavelength array
    """
    x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
    wavelength = slit.meta.wcs(x, y)[2]
    return wavelength


def radec_to_source_ids(catalog, source_ids=None, source_ra=None, source_dec=None, max_sep=1.0):
    """
    Convert source RA/Dec lists to source IDs from the catalog.

    If a ``source_ids`` list is provided, it will be combined with the
    source IDs found from the RA/Dec lists to form a union.

    Parameters
    ----------
    catalog : str
        The filename of the source catalog.

    source_ids : list
        List of source IDs to extract.

    source_ra : list of float
        Source right ascensions to be processed. The nearest matching source to each RA/Dec pair
        will be extracted. If both ``source_ids`` and ``source_ra``/``source_dec`` are provided,
        the lists will be combined and their union extracted.

    source_dec : list of float
        Source declinations to be processed, must have same length as ``source_ra``.

    max_sep : float
        Maximum separation in arcsec to consider a catalog source a match to the provided RA/Dec.

    Returns
    -------
    source_ids : ndarray or None
        List of unique source IDs to extract.
    """
    catalog = read_source_catalog(catalog)
    catalog_coord = catalog["sky_centroid"]
    if source_ids is None:
        source_ids = []
    else:
        # force_list coming into the step makes these all strings
        source_ids = np.atleast_1d(source_ids).astype(int).tolist()

    # check validity of RA/Dec inputs
    if source_ra is None and source_dec is not None:
        raise ValueError("source_ra must be provided if source_dec is provided.")
    if source_dec is None and source_ra is not None:
        raise ValueError("source_dec must be provided if source_ra is provided.")
    if (source_ra is not None) and (source_dec is not None):
        # force_list coming into the step makes these all strings
        source_ra = np.atleast_1d(source_ra).astype(float)
        source_dec = np.atleast_1d(source_dec).astype(float)
        if len(source_ra) != len(source_dec):
            raise ValueError("source_ra and source_dec must have the same length.")

        # find nearest catalog source for each RA/Dec pair
        for ra, dec in zip(source_ra, source_dec, strict=True):
            this_coord = SkyCoord(ra=ra, dec=dec, unit="deg")
            idx, sep, _dist3d = this_coord.match_to_catalog_sky(catalog_coord)
            if sep.arcsecond > max_sep:
                log.warning(
                    f"No catalog source found within {max_sep} arcsec of RA: {ra}, Dec: {dec}."
                )
                continue
            src_id = catalog["label"][idx]
            source_ids.append(src_id)

    if source_ids:
        return np.unique(np.atleast_1d(source_ids))  # return unique IDs only
    if source_ra is not None or source_dec is not None:
        raise ValueError(
            "source_ra and source_dec were provided, but no sources were found "
            "within source_max_sep of the requested location."
        )
    return None
