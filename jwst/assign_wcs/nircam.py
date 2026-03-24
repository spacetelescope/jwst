import logging

import asdf
import gwcs.coordinate_frames as cf
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import bind_bounding_box
from astropy.modeling.models import Const1D, Identity, Mapping, Shift
from gwcs import selector
from stdatamodels.jwst.datamodels import (
    DistortionModel,
    ImageModel,
    NIRCAMGrismModel,
    RegionsModel,
)
from stdatamodels.jwst.transforms.models import (
    NIRCAMBackwardGrismDispersion,
    NIRCAMForwardColumnGrismDispersion,
    NIRCAMForwardRowGrismDispersion,
)

from jwst.assign_wcs import pointing
from jwst.assign_wcs.util import (
    bounding_box_from_subarray,
    not_implemented_mode,
    subarray_transform,
    substripe_subarray_transforms,
    transform_bbox_from_shape,
    velocity_correction,
)
from jwst.lib import reffile_utils

log = logging.getLogger(__name__)


__all__ = ["create_pipeline", "dhs", "imaging", "tsgrism", "wfss"]


def create_pipeline(input_model, reference_files):
    """
    Create the WCS pipeline based on EXP_TYPE.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.wcs.WCS`.
    """
    log.debug(f"reference files used in NIRCAM WCS pipeline: {reference_files}")
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def imaging(input_model, reference_files):
    """
    Create the WCS pipeline for NIRCAM imaging data.

    It includes three coordinate frames - "detector", "v2v3", and "world"

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.ImageModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and filteroffset' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.wcs.WCS`.
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
    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        distortion = subarray2full | distortion

        # Bind the bounding box to the distortion model using the bounding box ordering
        # used by GWCS. This makes it clear the bounding box is set correctly to GWCS
        bind_bounding_box(distortion, bounding_box_from_subarray(input_model, order="F"), order="F")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )

    tel2sky = pointing.v23tosky(input_model)
    pipeline = [(detector, distortion), (v2v3, va_corr), (v2v3vacorr, tel2sky), (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the "detector" to "v2v3" transform for imaging mode.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.ImageModel` or \
                  `~stdatamodels.jwst.datamodels.CubeModel`
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
    transform = dist.model

    try:
        # Purposefully grab the bounding box tuple from the transform model in the
        # GWCS ordering
        bbox = transform.bounding_box.bounding_box(order="F")
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
        row = reffile_utils.find_row(filters, match_keys)
        if row is not None:
            col_offset = row.get("column_offset", "N/A")
            row_offset = row.get("row_offset", "N/A")
            log.debug(f"Offsets from filteroffset file are {col_offset}, {row_offset}")
            if col_offset != "N/A" and row_offset != "N/A":
                transform = Shift(col_offset) & Shift(row_offset) | transform
        else:
            log.debug("No match in fitleroffset file.")

    # Bind the bounding box to the distortion model using the bounding box ordering used by GWCS.
    # This makes it clear the bounding box is set correctly to GWCS
    bind_bounding_box(
        transform,
        transform_bbox_from_shape(input_model.data.shape, order="F") if bbox is None else bbox,
        order="F",
    )

    return transform


def _load_grism_models(specwcs_path):
    """Load dispersion models from a NIRCAMGrismModel reference file."""  # numpydoc ignore=RT01
    with NIRCAMGrismModel(specwcs_path) as f:
        displ = f.displ.instance
        dispx = f.dispx.instance
        dispy = f.dispy.instance
        invdispx = f.invdispx.instance if f.invdispx is not None else None
        invdispy = f.invdispy.instance if f.invdispy is not None else None
        invdispl = f.invdispl.instance if f.invdispl is not None else None
        orders = f.orders.instance
    return displ, dispx, dispy, invdispl, invdispx, invdispy, orders


def _build_grism_det2det(
    orders,
    displ,
    dispx,
    dispy,
    inv_lmodels=None,
    inv_xmodels=None,
    inv_ymodels=None,
    pupil="GRISMR",
):
    """Build the forward and backward grism dispersion transform pair."""  # numpydoc ignore=RT01
    if pupil == "GRISMC":
        forward_cls = NIRCAMForwardColumnGrismDispersion
    else:
        forward_cls = NIRCAMForwardRowGrismDispersion
    det2det = forward_cls(
        orders,
        lmodels=displ,
        xmodels=dispx,
        ymodels=dispy,
        inv_lmodels=inv_lmodels,
        inv_xmodels=inv_xmodels,
        inv_ymodels=inv_ymodels,
    )
    det2det.inverse = NIRCAMBackwardGrismDispersion(
        orders,
        lmodels=displ,
        xmodels=dispx,
        ymodels=dispy,
        inv_lmodels=inv_lmodels,
        inv_xmodels=inv_xmodels,
        inv_ymodels=inv_ymodels,
    )
    return det2det


def _apply_velocity_correction(det2det, velosys):
    """Apply barycentric velocity correction to the dispersion transform."""  # numpydoc ignore=RT01
    if velosys is not None:
        velocity_corr = velocity_correction(velosys)
        log.info(f"Added Barycentric velocity correction: {velocity_corr[1].amplitude.value}")
        det2det = det2det | Mapping((0, 1, 2, 3)) | Identity(2) & velocity_corr & Identity(1)
    return det2det


def _build_sky_pipeline_steps(input_model, reference_files, n_passthrough):
    """
    Build the distortion, DVA, and sky coordinate transforms for DHS and tsgrism modes.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'filteroffset' reference files.
    n_passthrough : int
        Number of extra dimensions to pass through unchanged alongside v2/v3.
        Use 2 for tsgrism (lam, order) and 3 for dhs (lam, order, stripe).

    Returns
    -------
    distortion, va_corr, tel2sky : astropy models
        The three transforms ready to slot into the WCS pipeline.
    """
    setra = Const1D(input_model.meta.wcsinfo.ra_ref)
    setra.inverse = Const1D(input_model.meta.wcsinfo.ra_ref)
    setdec = Const1D(input_model.meta.wcsinfo.dec_ref)
    setdec.inverse = Const1D(input_model.meta.wcsinfo.dec_ref)

    distortion = imaging_distortion(input_model, reference_files) & Identity(n_passthrough)

    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    ) & Identity(n_passthrough)

    tel2sky = pointing.v23tosky(input_model) & Identity(n_passthrough)

    return distortion, va_corr, tel2sky


def tsgrism(input_model, reference_files):
    """
    Create WCS pipeline for a NIRCAM Time Series Grism observation.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.CubeModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion', 'filteroffset', and 'specwcs' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.wcs.WCS`.

    Notes
    -----
    The TSGRISM mode should function effectively like the grism mode
    except that subarrays will be allowed. Since the transform models
    depend on the original full frame coordinates of the observation,
    the regular grism transforms will need to be shifted to the full
    frame coordinates around the trace transform.

    TSGRISM is only slated to work with GRISMR or DHS pupil elements and Module A.

    For this mode, the source is typically at crpix1 x crpix2, which
    are stored in keywords XREF_SCI, YREF_SCI.
    offset special requirements may be encoded in the X_OFFSET parameter,
    but those are handled in extract_2d.
    """
    # make sure this is a grism image
    if "NRC_TSGRISM" != input_model.meta.exposure.type:
        raise ValueError("The input exposure is not a NIRCAM time series grism")

    if input_model.meta.instrument.module != "A":
        raise ValueError("NRC_TSGRISM mode only supports module A")

    if "DHS" in input_model.meta.subarray.name:
        log.info("Building WCS for DHS data.")
        return dhs(input_model, reference_files)
    elif "GRISMR" not in input_model.meta.instrument.pupil:
        raise ValueError("NRC_TSGRISM mode only supports GRISMR and DHS pupil values.")

    frames = create_coord_frames()

    # Load disperser parameters and build the forward/backward transform pair.
    displ, dispx, dispy, invdispl, invdispx, _invdispy, orders = _load_grism_models(
        reference_files["specwcs"]
    )
    det2det = _build_grism_det2det(
        orders, displ, dispx, dispy, inv_lmodels=invdispl, inv_xmodels=invdispx
    )
    det2det = _apply_velocity_correction(det2det, input_model.meta.wcsinfo.velosys)

    # input into the forward transform is x,y,x0,y0,order
    # where x,y is the pixel location in the grism image
    # and x0,y0 is the source location in the "direct" image.
    # For this mode (tsgrism), it is assumed that the source is
    # at the nominal aperture reference point, i.e.,
    # crpix1 <--> xref_sci and crpix2 <--> yref_sci
    # offsets in X are handled in extract_2d, e.g. if an offset
    # special requirement was specified in the APT.
    xc, yc = (input_model.meta.wcsinfo.siaf_xref_sci, input_model.meta.wcsinfo.siaf_yref_sci)

    if xc is None:
        raise ValueError("XREF_SCI is missing.")

    if yc is None:
        raise ValueError("YREF_SCI is missing.")

    xcenter = Const1D(xc)
    xcenter.inverse = Const1D(xc)
    ycenter = Const1D(yc)
    ycenter.inverse = Const1D(yc)

    # x, y, order in goes to transform to full array location and order
    # get the shift to full frame coordinates
    sub_trans = subarray_transform(input_model)
    if sub_trans is not None:
        sub2direct = (
            sub_trans & Identity(1)
            | Mapping((0, 1, 0, 1, 2))
            | (Identity(2) & xcenter & ycenter & Identity(1))
            | det2det
        )
    else:
        sub2direct = (
            Mapping((0, 1, 0, 1, 2)) | (Identity(2) & xcenter & ycenter & Identity(1)) | det2det
        )

    distortion, va_corr, tel2sky = _build_sky_pipeline_steps(
        input_model, reference_files, n_passthrough=2
    )

    pipeline = [
        (frames["grism_detector"], sub2direct),
        (frames["direct_image"], distortion),
        (frames["v2v3"], va_corr),
        (frames["v2v3vacorr"], tel2sky),
        (frames["world"], None),
    ]

    return pipeline


def dhs(input_model, reference_files):
    """
    Create WCS pipeline for a NIRCAM Time Series DHS Grism observation.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.CubeModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion', 'filteroffset' and 'specwcs' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.wcs.WCS`.
    """
    if reference_files["regions"] == "":
        raise FileNotFoundError("No regions reference file provided.")

    frames = create_coord_frames()

    with RegionsModel(reference_files["regions"]) as regs_model:
        # Get the shift to full frame coordinates from stripe coords
        # Used in final transform from initial input x, y, order
        sub_trans_dict = substripe_subarray_transforms(input_model, regs_model)

        if regs_model.regions.shape == input_model.data.shape[-2:]:
            regions = regs_model.regions
        else:
            sub_regs_model = reffile_utils.get_subarray_model(input_model, regs_model)
            regions = sub_regs_model.regions
            sub_regs_model.close()

    label_mapper = selector.LabelMapperArray(
        regions,
        inputs_mapping=Mapping(mapping=(0, 1), n_inputs=3),
        inputs=("x", "y", "order"),
    )
    label_mapper.inverse = selector.LabelMapper(
        inputs=("x0", "y0", "lam", "order", "stripe"),
        mapper=Identity(1),
        inputs_mapping=Mapping((4,)),
    )

    with NIRCAMGrismModel(reference_files["specwcs"]) as f:
        displ = f.displ.instance
        dispx = f.dispx.instance
        dispy = f.dispy.instance
        orders = f.orders.instance
        stripes = f.stripes.instance
        fieldpoints = f.fieldpoints.instance

    if "LONG" in input_model.meta.instrument.detector.upper():
        longflag = True
        # Because nrcalong DHS uses existing transforms, we need to mock
        # the new structure to allow both to pass through this method.
        subarray_stripenum = int(input_model.meta.subarray.name.split("STRIPE")[1][0])
        stripes = np.array(range(subarray_stripenum)) + 1
        displ = [displ] * subarray_stripenum
        dispx = [dispx] * subarray_stripenum
        dispy = [dispy] * subarray_stripenum
    else:
        longflag = False
        # Short-wavelength GrismModels contain stripe transforms for both fieldpoints.
        # Down-select the transform lists to the relevant entries.
        fp_mask = [f in input_model.meta.aperture.pps_name for f in fieldpoints]
        displ = [b for (a, b) in zip(fp_mask, displ, strict=True) if a]
        dispx = [b for (a, b) in zip(fp_mask, dispx, strict=True) if a]
        dispy = [b for (a, b) in zip(fp_mask, dispy, strict=True) if a]
        stripes = [b for (a, b) in zip(fp_mask, stripes, strict=True) if a]

    # Initialize transforms dictionary to store stripe IDs as keys, transform models as values.
    # Used in RegionsSelector to choose correct transform given a stripe ID.
    transforms = {}

    velosys = input_model.meta.wcsinfo.velosys
    for i, stripe in enumerate(stripes):
        det2det = _build_grism_det2det(orders, displ[i], dispx[i], dispy[i])
        det2det = _apply_velocity_correction(det2det, velosys)

        xc, yc = (
            input_model.meta.wcsinfo.siaf_xref_sci - 1,
            input_model.meta.wcsinfo.siaf_yref_sci - 1,
        )
        if xc is None:
            raise ValueError("XREF_SCI is missing.")
        if yc is None:
            raise ValueError("YREF_SCI is missing.")

        if longflag:
            xform_refx = xform_refx.inverse = xcenter = xcenter.inverse = Const1D(xc)
            xform_refy = xform_refy.inverse = ycenter = ycenter.inverse = Const1D(yc)

        else:
            xcenter = Const1D(xc)
            xcenter.inverse = Const1D(0.0)
            ycenter = Const1D(yc)
            ycenter.inverse = Const1D(0.0)

            xform_refx = xform_refy = xform_refx.inverse = xform_refy.inverse = Const1D(0.0)

        stripe_model = Const1D(stripe)
        stripe_model.inverse = Const1D(stripe)

        if sub_trans_dict[stripe] is not None:
            sub2direct = (
                sub_trans_dict[stripe] & Identity(1)
                | Mapping((0, 1, 0, 1, 2, 2))
                | (Identity(2) & xform_refx & xform_refy & Identity(2))
                | det2det & stripe_model
                | (xcenter & ycenter & Identity(3))
            )
        else:
            sub2direct = (
                Mapping((0, 1, 0, 1, 2, 2))
                | (Identity(2) & xform_refx & xform_refy & Identity(2))
                | det2det & stripe_model
                | (xcenter & ycenter & Identity(3))
            )

        transforms[stripe] = sub2direct

    stripe2det = selector.RegionsSelector(
        inputs=["x", "y", "order"],
        outputs=["x0", "y0", "lam", "order", "stripe"],
        label_mapper=label_mapper,
        selector=transforms,
    )

    distortion, va_corr, tel2sky = _build_sky_pipeline_steps(
        input_model, reference_files, n_passthrough=3
    )

    pipeline = [
        (frames["grism_detector"], stripe2det),
        (frames["direct_image"], distortion),
        (frames["v2v3"], va_corr),
        (frames["v2v3vacorr"], tel2sky),
        (frames["world"], None),
    ]

    return pipeline


def wfss(input_model, reference_files):
    """
    Create the WCS pipeline for a NIRCAM grism observation.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.ImageModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion', 'filteroffset', and 'specwcs' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.wcs.WCS`.

    Notes
    -----
    The tree in the grism reference file has a section for each order.
    Not sure if there will be a separate passband reference file needed for
    the wavelength scaling or wedge offsets.

    The direct image the catalog has been created from was processed through
    resample, but the dispersed images have not been resampled. This is OK if
    the trace and dispersion solutions are defined with respect to the
    distortion-corrected image. The catalog from the combined direct image
    has object locations in in detector space and the RA DEC of the object on
    the sky.

    The WCS information for the grism image  plus the observed filter will be
    used to translate these to pixel locations for each of the objects.
    The grism images will then use their grism trace information to translate
    to detector space. The translation is assumed to be one-to-one for purposes
    of identifying the center of the object trace.

    The extent of the trace for each object can then be calculated based on
    the grism in use (row or column). Where the left/bottom of the trace starts
    at t = 0 and the right/top of the trace ends at t = 1, as long as they
    have been defined as such by the team.

    The extraction box is calculated to be the minimum bounding box of the
    object extent in the segmentation map associated with the direct image.
    The values of the min and max corners, taken from the computed minimum
    bounding box are saved in the photometry catalog in units of RA, DEC
    so they can be translated to pixels by the dispersed image's imaging-wcs.
    """
    # The input is the grism image
    if not isinstance(input_model, ImageModel):
        raise TypeError("The input data model must be an ImageModel.")

    # make sure this is a grism image
    if "NRC_WFSS" not in input_model.meta.exposure.type:
        raise ValueError("The input exposure is not a NIRCAM grism")

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
    displ, dispx, dispy, invdispl, invdispx, invdispy, orders = _load_grism_models(
        reference_files["specwcs"]
    )
    det2det = _build_grism_det2det(
        orders,
        displ,
        dispx,
        dispy,
        inv_lmodels=invdispl,
        inv_xmodels=invdispx,
        inv_ymodels=invdispy,
        pupil=input_model.meta.instrument.pupil,
    )

    # Add in the wavelength shift from the velocity dispersion
    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        velosys = None
    det2det = _apply_velocity_correction(det2det, velosys)

    # create the pipeline to construct a WCS object for the whole image
    # which can translate ra,dec to image frame reference pixels
    # it also needs to be part of the grism image wcs pipeline to
    # go from detector to world coordinates. However, the grism image
    # will be effectively translating pixel->world coordinates in a
    # manner that gives you the originating 'imaging' pixels ra and dec,
    # not the ra/dec on the sky from the pointing wcs of the grism image.
    image_pipeline = imaging(input_model, reference_files)

    # input is (x,y,x0,y0,order) -> x0, y0, wave, order
    # x0, y0 is in the image-detector reference frame already
    # and are fed to the wcs to calculate the ra,dec, pix offsets
    # and order are used to calculate the wavelength of the pixel
    grism_pipeline = [(gdetector, det2det)]

    # pass the x0,y0, wave, order, through the pipeline
    imagepipe = []
    world = image_pipeline.pop()[0]
    world.name = "sky"
    for cframe, trans in image_pipeline:
        trans = trans & (Identity(2))
        name = cframe.name
        cframe.name = name + "spatial"
        spatial_and_spectral = cf.CompositeFrame([cframe, spec], name=name)
        imagepipe.append((spatial_and_spectral, trans))

    # Output frame is Celestial + Spectral
    imagepipe.append((cf.CompositeFrame([world, spec], name="world"), None))
    grism_pipeline.extend(imagepipe)
    return grism_pipeline


def create_coord_frames():
    """
    Create the coordinate frames for NIRCAM imaging and grism modes.

    Returns
    -------
    frames : dict
        Dictionary of the coordinate frames.
    """
    gdetector = cf.Frame2D(name="grism_detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    detector = cf.Frame2D(
        name="full_detector", axes_order=(0, 1), axes_names=("dx", "dy"), unit=(u.pix, u.pix)
    )
    v2v3_spatial = cf.Frame2D(
        name="v2v3_spatial", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.deg, u.deg)
    )
    v2v3vacorr_spatial = cf.Frame2D(
        name="v2v3vacorr_spatial",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    sky_frame = cf.CelestialFrame(reference_frame=coord.ICRS(), name="icrs")
    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )
    frames = {
        "grism_detector": gdetector,
        "direct_image": cf.CompositeFrame([detector, spec], name="direct_image"),
        "v2v3": cf.CompositeFrame([v2v3_spatial, spec], name="v2v3"),
        "v2v3vacorr": cf.CompositeFrame([v2v3vacorr_spatial, spec], name="v2v3vacorr"),
        "world": cf.CompositeFrame([sky_frame, spec], name="world"),
    }
    return frames


exp_type2transform = {
    "nrc_image": imaging,
    "nrc_wfss": wfss,
    "nrc_tacq": imaging,
    "nrc_taconfirm": imaging,
    "nrc_coron": imaging,
    "nrc_focus": imaging,
    "nrc_tsimage": imaging,
    "nrc_tsgrism": tsgrism,
    "nrc_led": not_implemented_mode,
    "nrc_dark": not_implemented_mode,
    "nrc_flat": not_implemented_mode,
    "nrc_grism": not_implemented_mode,
}
