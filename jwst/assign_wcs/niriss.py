import logging
import numpy as np

import asdf
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import bind_bounding_box
from astropy.modeling.models import Const1D, Mapping, Identity, Shift
from astropy.modeling.bounding_box import CompoundBoundingBox
import gwcs.coordinate_frames as cf
from gwcs import wcs

from stdatamodels.jwst.datamodels import ImageModel, NIRISSGrismModel, DistortionModel
from stdatamodels.jwst.transforms.models import (
    NirissSOSSModel,
    NIRISSForwardRowGrismDispersion,
    NIRISSBackwardGrismDispersion,
    NIRISSForwardColumnGrismDispersion,
)

from .util import (
    not_implemented_mode,
    subarray_transform,
    velocity_correction,
    bounding_box_from_subarray,
    transform_bbox_from_shape,
)
from . import pointing
from jwst.lib.reffile_utils import find_row

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["create_pipeline", "imaging", "niriss_soss", "niriss_soss_set_input", "wfss"]


def create_pipeline(input_model, reference_files):
    """
    Create the WCS pipeline based on EXP_TYPE.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    log.debug(f"reference files used in NIRISS WCS pipeline: {reference_files}")
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def niriss_soss_set_input(model, order_number):
    """
    Extract a WCS for a specific spectral order.

    Parameters
    ----------
    model : `~jwst.datamodels.ImageModel`
        An instance of an ImageModel
    order_number : int
        The spectral order

    Returns
    -------
    gwcs.WCS
        The WCS corresponding to the spectral order.
    """
    # Make sure the spectral order is available.
    if order_number < 1 or order_number > 3:
        raise ValueError("Order must be between 1 and 3")

    # Return the correct transform based on the order_number
    obj = model.meta.wcs.forward_transform.get_model(order_number)

    # use the size of the input subarray
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )
    sky = cf.CelestialFrame(
        reference_frame=coord.ICRS(),
        axes_names=("ra", "dec"),
        axes_order=(0, 1),
        unit=(u.deg, u.deg),
        name="sky",
    )
    world = cf.CompositeFrame([sky, spec], name="world")
    pipeline = [(detector, obj), (world, None)]

    return wcs.WCS(pipeline)


def _niriss_order_bounding_box(input_model, order):
    bbox_y = np.array([-0.5, input_model.meta.subarray.ysize - 0.5])
    bbox_x = np.array([-0.5, input_model.meta.subarray.xsize - 0.5])

    if order == 1:
        return tuple(bbox_y), tuple(bbox_x)
    elif order == 2:
        return tuple(bbox_y), tuple(bbox_x)
    elif order == 3:
        return tuple(bbox_y), tuple(bbox_x)
    else:
        raise ValueError(
            f"Invalid spectral order: {order} provided. Spectral order must be 1, 2, or 3."
        )


def niriss_bounding_box(input_model):
    """
    Create a bounding box for the NIRISS model.

    Parameters
    ----------
    input_model : JwstDataModel
        The input datamodel.

    Returns
    -------
    CompoundBoundingBox
        The bounding box for the NIRISS model.
    """
    bbox = {(order,): _niriss_order_bounding_box(input_model, order) for order in [1, 2, 3]}
    model = input_model.meta.wcs.forward_transform
    return CompoundBoundingBox.validate(
        model, bbox, slice_args=[("spectral_order", True)], order="F"
    )


def niriss_soss(input_model, reference_files):
    """
    Create the WCS pipeline for NIRISS SOSS data.

    It includes TWO coordinate frames - "detector" and "world".

    Parameters
    ----------
    input_model : ImageModel, IFUImageModel, or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'specwcs' reference file.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    # Get the target RA and DEC, they will be used for setting the WCS RA
    # and DEC based on a conversation with Kevin Volk.
    try:
        target_ra = float(input_model.meta.target.ra)
        target_dec = float(input_model.meta.target.dec)
    except TypeError:
        # There was an error getting the target RA and DEC, so we are not going to continue.
        raise ValueError(
            f"Problem getting the TARG_RA or TARG_DEC from input model {input_model}"
        ) from None

    # Define the frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )
    sky = cf.CelestialFrame(
        reference_frame=coord.ICRS(),
        axes_names=("ra", "dec"),
        axes_order=(0, 1),
        unit=(u.deg, u.deg),
        name="sky",
    )
    world = cf.CompositeFrame([sky, spec], name="world")

    try:
        bytesio_or_path = reference_files["specwcs"]

        with asdf.open(bytesio_or_path) as af:
            wl1 = af.tree[1].copy()
            wl2 = af.tree[2].copy()
            wl3 = af.tree[3].copy()
    except Exception as e:
        raise OSError(
            f"Error reading wavelength correction from {reference_files['specwcs']}"
        ) from e

    velosys = input_model.meta.wcsinfo.velosys
    if velosys is not None:
        velocity_corr = velocity_correction(velosys)
        wl1 = wl1 | velocity_corr
        wl2 = wl2 | velocity_corr
        wl2 = wl3 | velocity_corr
        log.info(f"Applied Barycentric velocity correction: {velocity_corr[1].amplitude.value}")

    # Reverse the order of inputs passed to Tabular because it's in python order in modeling.
    # Consider changing it in modeling ?
    # fmt: off
    cm_order1 = (Mapping((0, 1, 1, 0)) |
        (Const1D(target_ra) & Const1D(target_dec) & wl1)
        ).rename("Order1")
    cm_order2 = (Mapping((0, 1, 1, 0)) |
        (Const1D(target_ra) & Const1D(target_dec) & wl2)
        ).rename("Order2")
    cm_order3 = (Mapping((0, 1, 1, 0)) |
        (Const1D(target_ra) & Const1D(target_dec) & wl3)
        ).rename("Order3")
    # fmt: on

    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        cm_order1 = subarray2full | cm_order1
        cm_order2 = subarray2full | cm_order2
        cm_order3 = subarray2full | cm_order3

        bbox = (
            (-0.5, input_model.meta.subarray.xsize - 0.5),
            (-0.5, input_model.meta.subarray.ysize - 0.5),
        )
        bind_bounding_box(cm_order1, bbox, order="F")
        bind_bounding_box(cm_order2, bbox, order="F")
        bind_bounding_box(cm_order3, bbox, order="F")

    # Define the transforms, they should accept (x,y) and return (ra, dec, lambda)
    soss_model = NirissSOSSModel([1, 2, 3], [cm_order1, cm_order2, cm_order3]).rename(
        "3-order SOSS Model"
    )

    # Define the pipeline based on the frames and models above.
    pipeline = [(detector, soss_model), (world, None)]

    return pipeline


def imaging(input_model, reference_files):
    """
    Create the WCS pipeline for NIRISS imaging data.

    It includes three coordinate frames - "detector" "v2v3" and "world".

    Parameters
    ----------
    input_model : ImageModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' reference file.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.arcsec, u.arcsec)
    )
    v2v3vacorr = cf.Frame2D(
        name="v2v3vacorr", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.arcsec, u.arcsec)
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    distortion = imaging_distortion(input_model, reference_files)

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )

    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        distortion = subarray2full | distortion

        # Bind the bounding box to the distortion model using the bounding box ordering
        # used by GWCS. This makes it clear the bounding box is set correctly to GWCS
        bind_bounding_box(distortion, bounding_box_from_subarray(input_model, order="F"), order="F")

    tel2sky = pointing.v23tosky(input_model)
    pipeline = [(detector, distortion), (v2v3, va_corr), (v2v3vacorr, tel2sky), (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the transform from "detector" to "v2v3".

    Parameters
    ----------
    input_model : ImageModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and filteroffset' reference files.

    Returns
    -------
    distortion : `astropy.modeling.Model`
        The transform from "detector" to "v2v3".
    """
    dist = DistortionModel(reference_files["distortion"])
    distortion = dist.model

    try:
        bbox = distortion.bounding_box.bounding_box(order="F")
    except NotImplementedError:
        # Check if the transform in the reference file has a ``bounding_box``.
        # If not set a ``bounding_box`` equal to the size of the image after
        # assembling all distortion corrections.
        bbox = None
    dist.close()

    # Add an offset for the filter
    if reference_files["filteroffset"] is not None:
        obsfilter = input_model.meta.instrument.filter
        obspupil = input_model.meta.instrument.pupil
        with asdf.open(reference_files["filteroffset"]) as filter_offset:
            filters = filter_offset.tree["filters"]

        match_keys = {"filter": obsfilter, "pupil": obspupil}
        row = find_row(filters, match_keys)
        if row is not None:
            col_offset = row.get("column_offset", "N/A")
            row_offset = row.get("row_offset", "N/A")
            log.info(f"Offsets from filteroffset file are {col_offset}, {row_offset}")

            if col_offset != "N/A" and row_offset != "N/A":
                distortion = Shift(col_offset) & Shift(row_offset) | distortion
        else:
            log.debug("No match in fitleroffset file.")

    # Bind the bounding box to the distortion model using the bounding box ordering used by GWCS.
    # This makes it clear the bounding box is set correctly to GWCS
    bind_bounding_box(
        distortion,
        transform_bbox_from_shape(input_model.data.shape, order="F") if bbox is None else bbox,
        order="F",
    )

    return distortion


def wfss(input_model, reference_files):
    """
    Create the WCS pipeline for a NIRISS grism observation.

    Parameters
    ----------
    input_model : ImageModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'specwcs' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.

    Notes
    -----
    The tree in the grism reference file has a section for each order/beam as
    well as the link to the filter data file, not sure if there will be a
    separate passband reference file needed for the wavelength scaling or the
    wedge offsets. This file is currently created in
    jwreftools/niriss/niriss_reftools.

    The direct image the catalog has been created from was corrected for
    distortion, but the dispersed images have not. This is OK if the trace and
    dispersion solutions are defined with respect to the distortion-corrected
    image. The catalog from the combined direct image has object locations in
    in detector space and the RA DEC of the object on sky.

    The WCS information for the grism image  plus the observed filter will be
    used to translate these to pixel locations for each of the objects.
    The grism images will then use their grism trace information to translate
    to detector space. The translation is assumed to be one-to-one for purposes
    of identifying the center of the object trace.

    The extent of the trace for each object can then be calculated based on
    the grism in use (row or column). Where the left/bottom of the trace starts
    at t = 0 and the right/top of the trace ends at t = 1, as long as they
    have been defined as such by th team.

    The extraction box is calculated to be the minimum bounding box of the
    object extent in the segmentation map associated with the direct image.
    The values of the min and max corners are saved in the photometry
    catalog in units of RA,DEC so they can be translated to pixels by
    the dispersed image's imaging wcs.

    The sensitivity information from the original aXe style configuration
    file needs to be modified by the passband of the filter used for
    the direct image to get the min and max wavelengths
    which correspond to t=0 and t=1, this currently has been done by the team
    and the min and max wavelengths to use to calculate t are stored in the
    grism reference file as wrange, which can be selected by wrange_selector
    which contains the filter names.

    Source catalog use moved to extract_2d.
    """
    # The input is the grism image
    if not isinstance(input_model, ImageModel):
        raise TypeError("The input data model must be an ImageModel.")

    # make sure this is a grism image
    if "NIS_WFSS" != input_model.meta.exposure.type:
        raise ValueError("The input exposure is not NIRISS grism")

    # Create the empty detector as a 2D coordinate frame in pixel units
    gdetector = cf.Frame2D(
        name="grism_detector",
        axes_order=(0, 1),
        axes_names=("x_grism", "y_grism"),
        unit=(u.pix, u.pix),
    )
    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )

    # translate the x,y detector-in to x,y detector out coordinates
    # Get the disperser parameters which are defined as a model for each
    # spectral order
    with NIRISSGrismModel(reference_files["specwcs"]) as f:
        dispx = f.dispx
        dispy = f.dispy
        displ = f.displ
        invdispl = f.invdispl
        orders = f.orders
        fwcpos_ref = f.fwcpos_ref

    # This is the actual rotation from the input model
    fwcpos = input_model.meta.instrument.filter_position
    if fwcpos is None:
        raise ValueError("FWCPOS keyword value not found in input image")

    # sep the row and column grism models
    # In "DMS" orientation (same parity as the sky), the GR150C spectra
    # are aligned more closely with the rows, and the GR150R spectra are
    # aligned more closely with the columns.
    if input_model.meta.instrument.filter.endswith("C"):
        det2det = NIRISSForwardRowGrismDispersion(
            orders, lmodels=displ, xmodels=dispx, ymodels=dispy, theta=fwcpos - fwcpos_ref
        )
    elif input_model.meta.instrument.filter.endswith("R"):
        det2det = NIRISSForwardColumnGrismDispersion(
            orders, lmodels=displ, xmodels=dispx, ymodels=dispy, theta=fwcpos - fwcpos_ref
        )
    else:
        raise ValueError(f"FILTER keyword {input_model.meta.instrument.filter} is not valid.")

    backward = NIRISSBackwardGrismDispersion(
        orders, lmodels=invdispl, xmodels=dispx, ymodels=dispy, theta=-(fwcpos - fwcpos_ref)
    )
    det2det.inverse = backward

    # Add in the wavelength shift from the velocity dispersion
    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    if velosys is not None:
        velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
        log.info(f"Added Barycentric velocity correction: {velocity_corr[1].amplitude.value}")
        det2det = det2det | Mapping((0, 1, 2, 3)) | Identity(2) & velocity_corr & Identity(1)

    # create the pipeline to construct a WCS object for the whole image
    # which can translate ra,dec to image frame reference pixels
    # it also needs to be part of the grism image wcs pipeline to
    # go from detector to world coordinates. However, the grism image
    # will be effectively translating pixel->world coordinates in a
    # manner that gives you the originating pixels ra and dec, not the
    # pure ra/dec on the sky from the pointing wcs.

    # use the imaging_distortion reference file here
    image_pipeline = imaging(input_model, reference_files)

    # forward input is (x,y,lam,order) -> x, y
    # backward input needs to be the same ra, dec, lam, order -> x, y
    grism_pipeline = [(gdetector, det2det)]

    # pass through the wave, beam  and theta in the pipeline
    # Theta is a constant for each grism exposure and is in the
    # meta information for the input_model, pass it to the model
    # so the user doesn't have to

    imagepipe = []
    world = image_pipeline.pop()[0]
    world.name = "sky"
    for cframe, trans in image_pipeline:
        trans = trans & (Identity(2))
        name = cframe.name
        cframe.name = name + "spatial"
        spatial_and_spectral = cf.CompositeFrame([cframe, spec], name=name)
        imagepipe.append((spatial_and_spectral, trans))
    imagepipe.append((cf.CompositeFrame([world, spec], name="world"), None))
    grism_pipeline.extend(imagepipe)

    return grism_pipeline


exp_type2transform = {
    "nis_image": imaging,
    "nis_wfss": wfss,
    "nis_soss": niriss_soss,
    "nis_ami": imaging,
    "nis_tacq": imaging,
    "nis_taconfirm": imaging,
    "nis_focus": imaging,
    "nis_dark": not_implemented_mode,
    "nis_lamp": not_implemented_mode,
}
